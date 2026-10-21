from datetime import date, datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from swisstrace import (assert_raw_crv, lookup_calibration, swisstrace_correct,
                       swisstrace_date, swisstrace_process, swisstrace_qc,
                       swisstrace_convert_batch)
from swisstrace.correction import HALF_LIVES


def test_exact_correction_and_metadata(raw_factory):
    path = raw_factory()
    res = swisstrace_correct(path, .5, "F18", pet_start=40)
    assert res["background"] == 30 and res["n_background"] == 40
    assert res["n_raw"] == 301 and res["acq_duration"] == 300
    assert res["t0_seconds"] == 40 and not res["t0_detected"]
    assert res["tac"].activity.iloc[0] == 0
    np.testing.assert_allclose(res["tac"].activity.iloc[1:], 50, atol=1e-12)
    np.testing.assert_allclose(res["tac"].frame_start.iloc[1:], res["tac"].frame_end.iloc[:-1])
    assert list(res["raw"]) == ["time", "coincidence", "singles1", "singles2"]
    assert list(res["tac"]) == ["time", "activity", "frame_start", "frame_end", "frame_dur"]
    assert swisstrace_date(path) == date(2026, 1, 1)


@pytest.mark.parametrize("pet_start", [40, "00:00:40", datetime(2026, 1, 1, 0, 0, 40, tzinfo=timezone.utc)])
def test_start_formats(raw_factory, pet_start):
    res = swisstrace_correct(raw_factory(), .5, half_life=6586.2, pet_start=pet_start,
                             zero_first_frame=False)
    np.testing.assert_allclose(res["tac"].activity, 50)
    assert res["isotope"] is None


@pytest.mark.parametrize("isotope,half_life", HALF_LIVES.items())
def test_isotopes(raw_factory, isotope, half_life):
    res = swisstrace_correct(raw_factory(), .5, isotope.lower(), pet_start=40)
    assert res["half_life"] == half_life


def test_auto_baseline(raw_factory):
    res = swisstrace_correct(raw_factory(), .5, "f-18")
    assert res["t0_seconds"] == 20
    assert res["lead"] == 20 and res["background"] == 30
    assert res["n_background"] == 40 and res["t0_detected"]


def test_auto_excludes_isolated_spike(raw_factory):
    counts = np.r_[np.full(40, 30.), np.full(30, 130.)]
    counts[20] = 1000
    res = swisstrace_correct(raw_factory(counts, np.arange(70)), .5, "F18")
    assert res["background"] == 30 and res["n_background"] == 39
    assert res["t0_seconds"] == 20


def test_short_flat_curve_and_early_onset(raw_factory):
    with pytest.warns(UserWarning, match="No sustained rise"):
        res = swisstrace_correct(raw_factory([30, 30], [0, 1]), .5, "F18")
    np.testing.assert_array_equal(res["tac"].activity, [0, 0])
    with pytest.warns(UserWarning):
        res = swisstrace_correct(raw_factory(np.r_[np.full(14, 30), np.full(10, 130)], np.arange(24)), .5, "F18")
    assert res["lead"] == 14 and res["t0_seconds"] == 0


def test_frame_boundaries_and_empty_frames(raw_factory):
    # Values after background subtraction/decay chosen to be 1, 3, 5, 7.
    times = np.array([0, 1, 2, 4, 6.])
    counts = 30 + np.array([0, 1, 3, 5, 7]) * np.exp(-np.log(2) / 6586.2 * (times - 1))
    res = swisstrace_correct(raw_factory(counts, times), 1, "F18", pet_start=1,
                             zero_first_frame=False, frame_scheme={"width": [1, 2], "end": [3, 7]})
    np.testing.assert_allclose(res["tac"].activity, [1, 3, np.nan, 5, 7], equal_nan=True)
    np.testing.assert_allclose(res["tac"].time, [.5, 1.5, 2.5, 4, 6])


@pytest.mark.parametrize("kwargs,match", [({}, "isotope is required"), ({"isotope": "X99"}, "Unknown isotope"),
    ({"half_life": 0}, "half_life"), ({"isotope": "F18", "pet_start": "25:00"}, "clock time"),
    ({"isotope": "F18", "pet_start": 9999}, "after the end"),
    ({"isotope": "F18", "min_run": 1.5}, "integer"),
    ({"isotope": "F18", "frame_scheme": {"width": [0], "end": [5]}}, "frame width")])
def test_invalid_inputs(raw_factory, kwargs, match):
    with pytest.raises(ValueError, match=match):
        swisstrace_correct(raw_factory(), .5, **kwargs)


def test_processed_rejected(tmp_path):
    path = tmp_path / "corrected.crv"
    path.write_text("Corrected_&_calibrated_[kBq/cc]\n0.5 0.0\n")
    with pytest.raises(ValueError, match="corrected/processed"):
        assert_raw_crv(path)


def test_calibration_matching(raw_factory):
    dates, factors = ["2026-01-03", "2026-01-01", "2026-01-01"], [.6, .4, .5]
    assert lookup_calibration(dates, factors, date="2026/01/02").calibration_factor.iloc[0] == .5
    assert lookup_calibration(dates, factors, date="2026-01-02", method="nearest").calibration_factor.iloc[0] == .4
    assert lookup_calibration(dates, factors, filename=raw_factory(), method="exact").calibration_factor.iloc[0] == .5
    with pytest.warns(UserWarning):
        result = lookup_calibration(dates, factors, date=["2025-01-01", "2026-02-01"], max_gap=2)
    assert result.calibration_factor.isna().all()
    assert result.gap_days.iloc[1] == 29
    with pytest.raises(ValueError, match="exactly one"):
        lookup_calibration(dates, factors)
    with pytest.warns(UserWarning):
        assert lookup_calibration([], [], date="2026-01-01").calibration_factor.isna().all()
    assert lookup_calibration(dates, factors, date=swisstrace_correct(raw_factory(), .5, "F18")).calibration_factor.iloc[0] == .5


def test_exports_and_overwrite(raw_factory, tmp_path):
    path = raw_factory()
    paths = swisstrace_process(path, .5, "F18", pet_start=40, sub="sub-01", ses="ses-02")
    assert all(p.exists() for p in paths.values())
    assert paths["plot_png"].read_bytes().startswith(b"\x89PNG")
    tsv = pd.read_csv(paths["bids_tsv"], sep="\t")
    np.testing.assert_allclose(tsv.whole_blood_radioactivity.iloc[1:], 50)
    meta = json.loads(paths["bids_json"].read_text())
    assert meta["WholeBloodAvail"] and not meta["PlasmaAvail"]
    assert meta["whole_blood_radioactivity"]["Units"] == "kBq/mL"
    before = {k: p.read_bytes() for k, p in paths.items()}
    swisstrace_process(path, .7, "F18", pet_start=40, sub="01", ses="02", overwrite=False)
    assert before == {k: p.read_bytes() for k, p in paths.items()}
    figure = swisstrace_qc(swisstrace_correct(path, .5, "F18"))
    assert len(figure.axes) == 2


def test_bids_image_matching(raw_factory, tmp_path):
    bids = tmp_path / "dataset"
    other = bids / "sub-10" / "pet"
    other.mkdir(parents=True)
    (other / "sub-10_pet.nii.gz").touch()
    output = swisstrace_process(raw_factory(), .5, "F18", sub="1", bids_dir=bids)
    assert output["bids_tsv"].parent == bids / "sub-1" / "pet"
    pet = bids / "sub-01" / "ses-02" / "pet"
    pet.mkdir(parents=True)
    (pet / "sub-01_ses-02_trc-FDG_run-1_pet.nii.gz").touch()
    output = swisstrace_process(raw_factory(), .5, "F18", sub="01", ses="02", bids_dir=bids)
    assert output["bids_tsv"] == pet / "sub-01_ses-02_trc-FDG_run-1_recording-autosampler_blood.tsv"
    with pytest.raises(ValueError, match="without `sub`"):
        swisstrace_process(raw_factory(), .5, "F18", ses="01")


@pytest.mark.parametrize("extension", ["csv", "tsv", "xlsx"])
def test_batch_formats_and_lookup(raw_factory, tmp_path, extension):
    if extension == "xlsx":
        pytest.importorskip("openpyxl")
    path = raw_factory()
    manifest = pd.DataFrame({" FILENAME ": [path.name, "missing.crv", path.name],
                             "ISOTOPE": ["F18"] * 3, "SUB": ["01", None, None],
                             "pet_start": ["40", "40", "00:00:40"],
                             "calibration_factor": [None, None, .4]})
    file = tmp_path / f"manifest.{extension}"
    if extension == "xlsx":
        manifest.to_excel(file, index=False)
    else:
        manifest.to_csv(file, sep="," if extension == "csv" else "\t", index=False)
    result = swisstrace_convert_batch(file, cal_dates=["2025-12-31"], cal_values=[.5])
    assert list(result.status) == ["ok", "error", "ok"]
    assert (tmp_path / "BIDS/sub-01/pet/sub-01_recording-autosampler_blood.tsv").exists()
    log = pd.read_csv(tmp_path / "calibration_factors.tsv", sep="\t")
    assert len(log) == 2 and log.gap_days.iloc[0] == 1


def test_batch_half_life_and_empty(raw_factory, tmp_path):
    result = swisstrace_convert_batch(pd.DataFrame({"filename": [str(raw_factory())],
                                      "half_life": [6586.2], "calibration_factor": [.5]}), output_folder=tmp_path)
    assert list(result.status) == ["ok"]
    result = swisstrace_convert_batch(pd.DataFrame(columns=["filename"]), output_folder=tmp_path)
    assert result.empty
