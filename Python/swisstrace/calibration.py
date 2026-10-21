"""Date-based calibration lookup."""
from datetime import date as Date, datetime
from itertools import islice
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd

from .correction import _timestamp


def swisstrace_date(file, tz="UTC"):
    """Read the acquisition date from the first timestamp within 50 lines."""
    with Path(file).open(encoding="utf-8") as stream:
        for line in islice(stream, 50):
            try:
                values = [float(v) for v in line.split()[:6]]
                if len(values) == 6 and np.isfinite(values).all():
                    return _timestamp(values, tz).date()
            except ValueError:
                continue
    raise ValueError(f"Could not read a timestamp from '{file}'.")


def _sequence(value):
    if isinstance(value, (str, Path, Date, datetime)) or value is None or np.isscalar(value):
        return [value]
    return list(value)


def _date(value):
    if value is None or pd.isna(value):
        return None
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, Date):
        return value
    if isinstance(value, str) and re.fullmatch(r"\d{4}[-/]\d{1,2}[-/]\d{1,2}", value.strip()):
        return Date(*(int(v) for v in re.split("[-/]", value.strip())))
    raise ValueError(f"Could not interpret '{value}' as a date; use YYYY-MM-DD or YYYY/MM/DD.")


def lookup_calibration(cal_dates, cal_values, filename=None, date=None,
                       method="last", max_gap=None, tz="UTC"):
    """Match each study to last/nearest/exact calibration; return a DataFrame.

    Supply exactly one of filename (one or many raw paths) or date (date values,
    date strings, or a correction result). Duplicate last/exact dates choose the
    last entry; nearest ties choose the first chronologically sorted entry.
    max_gap masks only the factor, preserving matched date and signed gap_days.
    """
    if (filename is None) == (date is None):
        raise ValueError("Supply exactly one of filename or date.")
    if method not in ("last", "nearest", "exact"):
        raise ValueError("method must be last, nearest, or exact.")
    if max_gap is not None and (not np.isfinite(max_gap) or max_gap < 0):
        raise ValueError("max_gap must be finite and nonnegative.")
    if isinstance(date, dict):
        date = date["date"]
    studies = ([swisstrace_date(f, tz) for f in _sequence(filename)] if filename is not None
               else [_date(d) for d in _sequence(date)])
    dates = [_date(d) for d in _sequence(cal_dates)]
    values = np.asarray(_sequence(cal_values), dtype=float)
    if len(dates) != len(values):
        raise ValueError("cal_dates and cal_values must have the same length.")
    if any(d is None for d in dates):
        raise ValueError("cal_dates contains unparseable dates.")
    order = sorted(range(len(dates)), key=lambda i: dates[i])
    dates, values = [dates[i] for i in order], values[order]
    rows = []
    for study in studies:
        candidates = [] if study is None else [i for i, d in enumerate(dates)
            if method == "nearest" or (method == "last" and d <= study) or (method == "exact" and d == study)]
        if not candidates:
            rows.append((study, np.nan, None, np.nan))
            continue
        index = (min(candidates, key=lambda i: abs((study - dates[i]).days))
                 if method == "nearest" else candidates[-1])
        gap = (study - dates[index]).days
        factor = values[index] if max_gap is None or abs(gap) <= max_gap else np.nan
        rows.append((study, factor, dates[index], gap))
    result = pd.DataFrame(rows, columns=["study_date", "calibration_factor", "cal_date", "gap_days"])
    missing = result.calibration_factor.isna().sum()
    if missing:
        warnings.warn(f"{missing} of {len(rows)} studies had no calibration factor within range.", stacklevel=2)
    return result
