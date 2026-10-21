"""Run both implementations on identical input when R + tibble/jsonlite exist."""
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from swisstrace import swisstrace_correct, swisstrace_process, lookup_calibration

R_SCRIPT = r'''
args <- commandArgs(trailingOnly = TRUE)
for (f in list.files(args[1], pattern = "[.]R$", full.names = TRUE)) source(f)
mode <- args[4]
a <- list(file = args[2], calibration_factor = 0.425, isotope = "F18")
if (mode != "auto") a$pet_start <- 40.25
if (mode == "framed") a$frame_scheme <- data.frame(width = c(1, 10), end = c(180, 600))
if (mode == "nozero") a$zero_first_frame <- FALSE
if (mode == "c11") a$isotope <- "C11"
res <- do.call(swisstrace_correct, a)
jsonlite::write_json(res[c("tac", "raw", "background", "t0_seconds", "n_background",
                         "half_life", "lambda", "lead", "n_raw")], args[3],
                    digits = NA, auto_unbox = TRUE, na = "null")
'''


@pytest.fixture(scope="module")
def rscript():
    executable = shutil.which("Rscript")
    if executable is None:
        pytest.skip("Rscript unavailable; R parity tests run in CI")
    check = subprocess.run([executable, "-e", 'stopifnot(requireNamespace("tibble", quietly=TRUE), requireNamespace("jsonlite", quietly=TRUE))'], capture_output=True)
    if check.returncode:
        pytest.skip("R tibble/jsonlite unavailable")
    return executable


@pytest.mark.parametrize("mode", ["native", "auto", "framed", "nozero", "c11"])
def test_r_numerical_parity(raw_factory, tmp_path, rscript, mode):
    times = np.arange(0, 301, .5)
    counts = 30 + 2 * np.sin(times)
    post = times >= 60
    counts[post] += 150 * (1 - np.exp(-(times[post] - 60) / 2)) * np.exp(-(times[post] - 60) / 300)
    path = raw_factory(counts, times)
    options = {} if mode == "auto" else dict(pet_start=40.25)
    if mode == "framed":
        options["frame_scheme"] = {"width": [1, 10], "end": [180, 600]}
    if mode == "nozero":
        options["zero_first_frame"] = False
    res = swisstrace_correct(path, .425, "C11" if mode == "c11" else "F18", **options)
    script = tmp_path / "parity.R"
    script.write_text(R_SCRIPT)
    output = tmp_path / "r.json"
    source = Path(__file__).resolve().parents[2] / "R" / "R"
    subprocess.run([rscript, str(script), str(source), str(path), str(output), mode], check=True, capture_output=True, text=True)
    reference = json.loads(output.read_text())
    for key in ("tac", "raw"):
        expected = pd.DataFrame(reference[key])
        assert list(expected) == list(res[key])
        np.testing.assert_allclose(res[key], expected, rtol=1e-10, atol=1e-10, equal_nan=True)
    for key in ("background", "t0_seconds", "n_background", "half_life", "lambda", "lead", "n_raw"):
        np.testing.assert_allclose(res[key], np.nan if reference[key] is None else reference[key], rtol=1e-10, atol=1e-10, equal_nan=True)


def test_r_export_and_lookup_parity(raw_factory, tmp_path, rscript):
    path = raw_factory()
    root = Path(__file__).resolve().parents[2]
    script = tmp_path / "exports.R"
    script.write_text(r'''
a <- commandArgs(trailingOnly=TRUE)
for (f in list.files(a[1], pattern="[.]R$", full.names=TRUE)) source(f)
swisstrace_process(a[2], 0.5, "F18", pet_start=40, output_folder=a[3], sub="01", ses="02")
x <- lookup_calibration(c("2025-12-31", "2026-01-02"), c(.4, .5), date=c("2026-01-01", "2026-01-03"))
write.table(x, file.path(a[3], "lookup.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
''')
    r_out = tmp_path / "R-output"
    subprocess.run([rscript, str(script), str(root / "R/R"), str(path), str(r_out)], check=True, capture_output=True, text=True)
    paths = swisstrace_process(path, .5, "F18", pet_start=40, output_folder=tmp_path / "python-output", sub="01", ses="02")
    for key, relative, skip in [("corrected_crv", "Corrected_PMOD/raw_corrected.crv", 1),
                               ("bids_tsv", "BIDS/sub-01/ses-02/pet/sub-01_ses-02_recording-autosampler_blood.tsv", 1)]:
        np.testing.assert_allclose(np.loadtxt(paths[key], skiprows=skip), np.loadtxt(r_out / relative, skiprows=skip), atol=1e-12)
    meta = json.loads(paths["bids_json"].read_text())
    expected = json.loads((r_out / "BIDS/sub-01/ses-02/pet/sub-01_ses-02_recording-autosampler_blood.json").read_text())
    for key in ("PlasmaAvail", "WholeBloodAvail", "MetaboliteAvail", "DispersionCorrected"):
        assert meta[key] == expected[key]
    actual = lookup_calibration(["2025-12-31", "2026-01-02"], [.4, .5], date=["2026-01-01", "2026-01-03"])
    reference = pd.read_csv(r_out / "lookup.tsv", sep="\t")
    np.testing.assert_allclose(actual[["calibration_factor", "gap_days"]], reference[["calibration_factor", "gap_days"]])
