"""Correction of raw twilite coincidence data, matching the R implementation."""
from datetime import datetime, timedelta
from itertools import islice
from pathlib import Path
import re
import warnings
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

HALF_LIVES = dict(F18=6586.2, C11=1223.4, N13=597.9, O15=122.24,
                 Ga68=4057.74, Cu62=584.4, Zr89=282240, I124=360806.4, Rb82=76.4)


def assert_raw_crv(file):
    """Reject empty or already corrected/processed recordings; return True."""
    with Path(file).open(encoding="utf-8") as stream:
        lines = [s.strip() for s in islice(stream, 5) if s.strip()]
    if not lines:
        raise ValueError(f"File '{file}' is empty.")
    try:
        numeric = np.isfinite(float(lines[0].split()[0]))
    except ValueError:
        numeric = False
    if not numeric or re.search(r"corrected|kBq/cc|value\[|time\[seconds\]", lines[0], re.I):
        raise ValueError(f"'{Path(file).name}' looks like a corrected/processed .crv, "
                         "not a raw twilite recording. Expected YYYY M D H M S "
                         "coincidences singles1 singles2.")
    return True


def _timestamp(values, tz):
    return datetime(*(int(v) for v in values[:5]), tzinfo=ZoneInfo(tz)) + timedelta(seconds=float(values[5]))


def _number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, (str, bool)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a single finite number.")
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a single finite number.") from None
    if not np.isfinite(value) or (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError(f"Invalid {name}: must be finite" + (" and positive." if positive else "."))
    return value


def swisstrace_correct(file, calibration_factor, isotope=None, pet_start=None,
                       half_life=None, lead=20, frame_scheme=None,
                       zero_first_frame=True, baseline_k=3, min_run=5,
                       baseline_init=10, baseline_min=20, tz="UTC"):
    """Return correction metadata and pandas ``tac``/``raw`` DataFrames in a dict.

    ``pet_start`` accepts seconds from acquisition start, a datetime, or HH:MM[:SS]
    on the acquisition date. None detects a sustained rise using median + 3 MAD
    (R's 1.4826 scale), then places time zero ``lead`` seconds before it.
    ``frame_scheme`` is a DataFrame or mapping with width/end columns (seconds).
    Native samples are retained by default. Framed values are sample means over
    left-closed, right-open intervals; empty frames contain NaN.
    """
    assert_raw_crv(file)
    calibration_factor = _number(calibration_factor, "calibration_factor", positive=True)
    if half_life is None:
        if not isinstance(isotope, str) or not isotope:
            raise ValueError("isotope is required unless half_life is supplied.")
        key = re.sub(r"[^A-Za-z0-9]", "", isotope).upper()
        lut = {k.upper(): v for k, v in HALF_LIVES.items()}
        if key not in lut:
            raise ValueError(f"Unknown isotope '{isotope}'; supply half_life or use {list(HALF_LIVES)}.")
        half_life = lut[key]
    half_life = _number(half_life, "half_life", positive=True)
    lead = _number(lead, "lead", nonnegative=True)
    baseline_k = _number(baseline_k, "baseline_k", nonnegative=True)
    baseline_init = _number(baseline_init, "baseline_init", positive=True)
    baseline_min = _number(baseline_min, "baseline_min", nonnegative=True)
    min_run = _number(min_run, "min_run", positive=True)
    if not min_run.is_integer():
        raise ValueError("min_run must be a positive integer.")
    min_run = int(min_run)
    rows = []
    with Path(file).open(encoding="utf-8") as stream:
        for line in stream:
            parts = line.split()
            if len(parts) >= 9:
                rows.append([float(v) for v in parts[:9]])
    if not rows:
        raise ValueError(f"No parseable data rows found in '{file}'.")
    data = np.asarray(rows)
    if not np.isfinite(data).all():
        raise ValueError("Raw timestamps and counts must be finite.")
    stamps = [_timestamp(row, tz) for row in data]
    secs = np.array([s.timestamp() - stamps[0].timestamp() for s in stamps])
    if np.any(np.diff(secs) < 0):
        raise ValueError("Raw timestamps must be in chronological order.")
    coinc = data[:, 6]
    n = len(coinc)
    lead_used = np.nan
    bg_idx = None
    if pet_start is not None:
        if isinstance(pet_start, str):
            clock = pet_start.strip()
            if not re.fullmatch(r"[0-9]{1,2}:[0-9]{2}(:[0-9]{2}(\.[0-9]+)?)?", clock):
                raise ValueError("pet_start is not a valid clock time; use HH:MM:SS.")
            parts = clock.split(":")
            hour, minute = map(int, parts[:2])
            second = float(parts[2]) if len(parts) == 3 else 0
            if hour > 23 or minute > 59 or second >= 60:
                raise ValueError("pet_start is not a valid clock time.")
            pet_start = stamps[0].replace(hour=hour, minute=minute, second=0, microsecond=0) + timedelta(seconds=second)
        if isinstance(pet_start, datetime):
            if pet_start.tzinfo is None:
                pet_start = pet_start.replace(tzinfo=ZoneInfo(tz))
            t0 = pet_start.timestamp() - stamps[0].timestamp()
        else:
            t0 = _number(pet_start, "pet_start")
    else:
        bg_idx = list(np.flatnonzero(secs < baseline_init))
        if len(bg_idx) < 3:
            bg_idx = list(range(min(3, n)))
        onset = None
        for i in range(max(bg_idx) + 1, n):
            base = coinc[bg_idx]
            median = np.median(base)
            spread = 1.4826 * np.median(np.abs(base - median))
            if spread == 0 or not np.isfinite(spread):
                spread = np.std(base, ddof=1) if len(base) > 1 else 0
            if spread == 0 or not np.isfinite(spread):
                spread = 1
            threshold = median + baseline_k * spread
            if i + min_run <= n and np.all(coinc[i:i + min_run] > threshold):
                onset = i
                break
            if coinc[i] <= threshold:
                bg_idx.append(i)
        if onset is None:
            warnings.warn("No sustained rise above background detected; using start of recording as time 0.", stacklevel=2)
            t0 = 0.
        else:
            onset_sec = secs[onset]
            if onset_sec < baseline_min:
                warnings.warn("Rise detected with less than baseline_min seconds; background may be unreliable.", stacklevel=2)
            lead_used = min(lead, onset_sec)
            if lead_used < lead:
                warnings.warn(f"Only {onset_sec:g} s before rise; using a {lead_used:g} s lead.", stacklevel=2)
            t0 = onset_sec - lead_used
    if t0 < 0:
        warnings.warn("Time 0 is before recording start; clamping to 0.", stacklevel=2)
        t0 = 0.
    if t0 > secs[-1]:
        raise ValueError("pet_start is after the end of the recording.")
    if bg_idx is None:
        bg_idx = np.flatnonzero(secs < t0)
    background = float(np.mean(coinc[bg_idx])) if len(bg_idx) else 0.
    if not len(bg_idx):
        warnings.warn("No samples before supplied time 0; background set to 0.", stacklevel=2)
    decay = np.log(2) / half_life
    tau = secs - t0
    keep = tau >= 0
    times = tau[keep]
    values = (coinc[keep] - background) * np.exp(decay * times) * calibration_factor
    if frame_scheme is None:
        if len(times) >= 2:
            mids = (times[1:] + times[:-1]) / 2
            starts = np.r_[0., mids]
            ends = np.r_[mids, 2 * times[-1] - mids[-1]]
        else:
            starts, ends = np.zeros(len(times)), times.copy()
    else:
        scheme = pd.DataFrame(frame_scheme)
        if scheme.empty or not {"width", "end"}.issubset(scheme.columns):
            raise ValueError("frame_scheme needs nonempty width and end columns.")
        edges, previous = [0.], 0.
        for width, end in scheme[["width", "end"]].itertuples(index=False, name=None):
            width = _number(width, "frame width", positive=True)
            end = _number(end, "frame end", positive=True)
            if end <= previous or end - previous < width:
                raise ValueError("Frame ends must increase by at least their width.")
            count = int(np.floor((end - previous) / width + 1e-12))
            edges.extend(previous + width * np.arange(1, count + 1))
            previous = end
        edges = np.unique(edges)
        nframes = min(int(np.searchsorted(edges, times[-1], side="right")), len(edges) - 1)
        bins = np.searchsorted(edges, times, side="right") - 1
        framed = np.full(nframes, np.nan)
        for index in range(nframes):
            selected = values[bins == index]
            if len(selected):
                framed[index] = np.mean(selected)
        starts, ends = edges[:nframes], edges[1:nframes + 1]
        times, values = (starts + ends) / 2, framed
    if zero_first_frame and len(values):
        values[0] = 0.
    return dict(tac=pd.DataFrame(dict(time=times, activity=values, frame_start=starts,
                                     frame_end=ends, frame_dur=ends - starts)),
                calibration_factor=calibration_factor, background=background,
                isotope=isotope, half_life=half_life, **{"lambda": decay},
                date=stamps[0].date(), start_time=stamps[0].strftime("%H:%M:%S.%f")[:-3],
                acq_start=stamps[0], acq_end=stamps[-1], acq_duration=float(secs[-1]),
                pet_start=stamps[0] + timedelta(seconds=float(t0)), t0_seconds=float(t0),
                t0_detected=pet_start is None, lead=lead_used, n_raw=n,
                n_background=len(bg_idx), file=str(Path(file).resolve()), frame_scheme=frame_scheme,
                raw=pd.DataFrame(dict(time=secs, coincidence=coinc, singles1=data[:, 7], singles2=data[:, 8])))
