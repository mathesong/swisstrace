"""Manifest-driven batch processing with per-row error reporting."""
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .calibration import lookup_calibration
from .process import swisstrace_process


def _read_manifest(manifest):
    if isinstance(manifest, pd.DataFrame):
        data, directory = manifest.copy(), Path.cwd()
    else:
        path = Path(manifest)
        directory = path.resolve().parent
        if path.suffix.lower() in (".xlsx", ".xls"):
            data = pd.read_excel(path, dtype=str)
        elif path.suffix.lower() in (".csv", ".tsv", ".tab", ".txt"):
            data = pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", dtype=str)
        else:
            raise ValueError("Unsupported manifest type; use .xlsx, .xls, .csv or .tsv.")
    data.columns = [str(c).strip().lower() for c in data.columns]
    if data.columns.duplicated().any():
        raise ValueError("Manifest has duplicate column names after case normalization.")
    if "filename" not in data:
        raise ValueError("Manifest is missing required column: filename.")
    return data, directory


def _cell(row, name):
    value = row.get(name)
    if value is None or pd.isna(value) or (isinstance(value, str) and not value.strip()):
        return None
    return value.strip() if isinstance(value, str) else value


def swisstrace_convert_batch(manifest, output_folder=None, cal_dates=None,
                            cal_values=None, bids_dir=None, cal_method="last",
                            cal_max_gap=None, **kwargs):
    """Process a CSV/TSV/Excel manifest or DataFrame and return row statuses.

    Columns: filename, isotope (or half_life), pet_start, calibration_factor,
    sub, ses. Relative filenames are resolved against the manifest directory.
    Missing factors use the calibration log and are recorded in
    calibration_factors.tsv. Per-row lookup and processing errors are isolated.
    Install the excel extra to read .xlsx/.xls. String labels preserve leading zeros.
    """
    if cal_method not in ("last", "nearest", "exact"):
        raise ValueError("cal_method must be last, nearest, or exact.")
    data, directory = _read_manifest(manifest)
    root = Path(output_folder) if output_folder is not None else directory
    factors = pd.to_numeric(data.get("calibration_factor", pd.Series(np.nan, index=data.index)), errors="coerce").to_numpy()
    if np.isnan(factors).any() and (cal_dates is None or cal_values is None):
        raise ValueError("Missing calibration_factor; supply cal_dates and cal_values.")
    root.mkdir(parents=True, exist_ok=True)
    statuses, log = [], []
    for i, (_, row) in enumerate(data.iterrows()):
        filename, isotope = _cell(row, "filename"), _cell(row, "isotope")
        factor, message, status = factors[i], None, "ok"
        logged = False
        try:
            if filename is None:
                raise ValueError("filename is missing.")
            path = Path(filename)
            if not path.is_absolute():
                path = directory / path
            if np.isnan(factor):
                match = lookup_calibration(cal_dates, cal_values, filename=path,
                                           method=cal_method, max_gap=cal_max_gap).iloc[0]
                factor = match.calibration_factor
                log.append(dict(filename=filename, date=match.study_date,
                                calibration_factor=factor, gap_days=match.gap_days))
                logged = True
            if not np.isfinite(factor):
                raise ValueError("No calibration factor found within range.")
            options = dict(kwargs)
            pet_start = _cell(row, "pet_start")
            if isinstance(pet_start, str) and re.fullmatch(r"[+-]?(\d+(\.\d*)?|\.\d+)", pet_start):
                pet_start = float(pet_start)
            options.update(pet_start=pet_start, sub=_cell(row, "sub"), ses=_cell(row, "ses"))
            half_life = _cell(row, "half_life")
            if half_life is not None:
                options["half_life"] = float(half_life)
            swisstrace_process(path, factor, isotope, output_folder=root,
                               bids_dir=bids_dir if options["sub"] is not None else None, **options)
        except Exception as exc:
            status, message = "error", str(exc)
            if np.isnan(factors[i]) and not logged:
                log.append(dict(filename=filename, date=None, calibration_factor=np.nan, gap_days=np.nan))
        statuses.append(dict(filename=filename, isotope=isotope, calibration_factor=factor,
                             status=status, message=message))
    if np.isnan(factors).any():
        pd.DataFrame(log, columns=["filename", "date", "calibration_factor", "gap_days"]).to_csv(
            root / "calibration_factors.tsv", sep="\t", index=False, na_rep="NA")
    return pd.DataFrame(statuses, columns=["filename", "isotope", "calibration_factor", "status", "message"])
