"""PMOD, QC and BIDS PET blood exports."""
import json
from pathlib import Path
import re
import warnings

import numpy as np

from .correction import swisstrace_correct
from .qc import swisstrace_qc


def _label(value, entity):
    label = re.sub(r"[^A-Za-z0-9]", "", re.sub(rf"^{entity}-", "", str(value)))
    if not label:
        raise ValueError(f"{entity} must contain an alphanumeric BIDS label.")
    return label


def _format(value):
    # Preserve full precision, including fractional acquisition times.
    return np.format_float_positional(float(value), unique=True, trim="0") if np.isfinite(value) else "n/a"


def swisstrace_process(filename, calibration_factor, isotope=None, *,
                       output_folder=None, sub=None, ses=None, bids_dir=None,
                       recording="autosampler", overwrite=True, **kwargs):
    """Correct a recording and return a dict of written pathlib Paths.

    Additional keyword arguments are forwarded to swisstrace_correct. With sub,
    write BIDS blood TSV/JSON; bids_dir locates a matching PET image and inherits
    its entities. Multiple matches warn and use the first sorted path, like R.
    overwrite=False preserves each existing file individually.
    """
    if sub is None and (ses is not None or bids_dir is not None):
        raise ValueError("ses or bids_dir was provided without `sub`.")
    if sub is not None:
        subject = _label(sub, "sub")
        session = _label(ses, "ses") if ses is not None else None
        recording = _label(recording, "recording")
    res = swisstrace_correct(filename, calibration_factor, isotope, **kwargs)
    root = Path(output_folder) if output_folder is not None else Path(filename).parent
    base = re.sub(r"\.crv$", "", Path(filename).name, flags=re.I)
    paths = dict(corrected_crv=root / "Corrected_PMOD" / f"{base}_corrected.crv",
                 plot_png=root / "Plots" / f"{base}_corrected.png")
    if sub is not None:
        bids = Path(bids_dir) if bids_dir is not None else root / "BIDS"
        directory = bids / f"sub-{subject}"
        entities = f"sub-{subject}"
        if session is not None:
            directory /= f"ses-{session}"
            entities += f"_ses-{session}"
        directory /= "pet"
        if bids_dir is not None:
            pattern = re.compile(rf"^{entities}(_.*)?_pet\.nii(\.gz)?$")
            hits = sorted(p for p in bids.rglob("*_pet.nii*") if p.is_file() and pattern.fullmatch(p.name))
            if hits:
                if len(hits) > 1:
                    warnings.warn(f"{len(hits)} PET images matched; using the first: {hits[0]}", stacklevel=2)
                directory = hits[0].parent
                entities = re.sub(r"_pet\.nii(\.gz)?$", "", hits[0].name)
        stem = f"{entities}_recording-{recording}_blood"
        paths.update(bids_tsv=directory / f"{stem}.tsv", bids_json=directory / f"{stem}.json")
    rows = "".join(f"{_format(t)}\t{_format(a)}\n" for t, a in
                   res["tac"][["time", "activity"]].itertuples(index=False, name=None))
    metadata = dict(PlasmaAvail=False, WholeBloodAvail=True, MetaboliteAvail=False,
                    DispersionCorrected=False,
                    time=dict(Description="Time relative to time zero (PET scan start); native sample time or frame midpoint.", Units="s"),
                    whole_blood_radioactivity=dict(Description="Whole blood radioactivity concentration: background-subtracted, decay-corrected to time zero, and calibrated.", Units="kBq/mL"))
    for kind, path in paths.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and not overwrite:
            continue
        if kind == "plot_png":
            figure = swisstrace_qc(res)
            figure.savefig(path, dpi=110)
            figure.clear()
        elif kind == "bids_json":
            path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        else:
            header = ("time\twhole_blood_radioactivity\n" if kind == "bids_tsv" else
                      "Corrected_&_calibrated_[kBq/cc]_>___time[seconds]\tvalue[kBq/cc]\n")
            path.write_text(header + rows, encoding="utf-8")
    return paths
