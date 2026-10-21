# swisstrace for Python

Python equivalents of the [R tools](../R/README.md), using NumPy, pandas and
Matplotlib. Python 3.10 or newer is required.

## Install

From the repository root:

```sh
python -m pip install ./Python
# Optional Excel manifests and development tests:
python -m pip install -e './Python[excel,test]'
```

## Correct and export

```python
from swisstrace import swisstrace_correct, swisstrace_process, swisstrace_qc

res = swisstrace_correct(
    "recording.crv", calibration_factor=0.425, isotope="F18",
    pet_start="11:16:30",
)
print(res["tac"])  # pandas: time, activity, frame_start, frame_end, frame_dur
fig = swisstrace_qc(res)
fig.savefig("qc.png")

paths = swisstrace_process(
    "recording.crv", calibration_factor=0.425, isotope="F18",
    pet_start="11:16:30", output_folder="converted", sub="01", ses="02",
)
```

Outputs are `Corrected_PMOD/<base>_corrected.crv`,
`Plots/<base>_corrected.png` and, with `sub`, BIDS blood TSV/JSON files under
`BIDS/sub-01/ses-02/pet/`. Set `bids_dir` to an existing dataset to write beside a
matching PET image, inheriting its entities. Multiple matching images warn and
use the first sorted path, matching the R behavior. `overwrite=False` preserves
each existing output. Session or dataset arguments require a subject.

## API correspondence

| R function | Python function / return |
| --- | --- |
| `swisstrace_correct()` | Same name; dictionary with `tac` and `raw` DataFrames plus metadata |
| `swisstrace_process()` | Same name; dictionary of `pathlib.Path` output paths |
| `swisstrace_qc()` | Same name; Matplotlib `Figure` (R returns its input invisibly) |
| `lookup_calibration()` | Same name; pandas DataFrame |
| `swisstrace_date()` | Same name; `datetime.date` |
| `swisstrace_convert_batch()` | Same name; pandas status DataFrame |
| `assert_raw_crv()` | Same name; `True`, or raises `ValueError` |

Use `None` for R's `NULL`, `True`/`False` for logical arguments, dictionary indexing
for R list elements, and DataFrames or column dictionaries for frame schemes.
`res["lambda"]` stores the decay constant. No R installation is needed to run Python.

## Correction options

`pet_start` accepts seconds after recording start, a Python datetime, or an
`HH:MM[:SS[.fraction]]` clock time on the acquisition date. A naive datetime is
interpreted using `tz` (default `UTC`). For a scan after midnight, supply a full
datetime or elapsed seconds. Omitting `pet_start` enables the same forward
median/MAD baseline detector as R, with `lead=20`, `baseline_k=3`, `min_run=5`,
`baseline_init=10` and `baseline_min=20` (times in seconds).

Known isotopes: F18, C11, N13, O15, Ga68, Cu62, Zr89, I124 and Rb82. Names are case
and punctuation insensitive. Supply `half_life` in seconds to override the table
or use an unlisted isotope. The coincidence channel alone is corrected; both
singles channels are preserved in `res["raw"]`.

```python
res = swisstrace_correct(
    "recording.crv", .425, "F18", pet_start=40,
    frame_scheme={"width": [1, 10], "end": [180, 600]},
    zero_first_frame=False,
)
```

With no scheme, the original sample times are retained. Framed values are means
of samples in left-closed, right-open intervals. Empty frames contain `NaN`;
missing exported values are `n/a`. Frames beyond acquired data are omitted.
`zero_first_frame=True` (default) forces the first sample/frame value to zero.

The Python port deliberately rejects nonfinite/nonpositive calibration factors
and half-lives, invalid frame schemes, out-of-order timestamps, and a scan start
after the recording ends. Short flat recordings return a warning and a baseline
estimate. These guards make invalid inputs explicit; valid-data mathematics and
defaults follow R. Native sample times are described accurately in the Python
BIDS sidecar; numerical values and availability flags match R.

## Calibration and batches

```python
from swisstrace import lookup_calibration, swisstrace_convert_batch

cal = lookup_calibration(
    ["2025-09-04", "2026-01-15", "2026-02-20"], [.198, .41, .425],
    filename="recording.crv",  # or date="2026-02-22", not both
    method="last", max_gap=90,
)
status = swisstrace_convert_batch(
    "manifest.csv", output_folder="converted",
    cal_dates=["2026-01-15", "2026-02-20"], cal_values=[.41, .425],
)
print(status[["filename", "status", "message"]])
```

Lookup supports scalar or multiple study filenames/dates and a correction result
as `date`. Methods: `last` (latest on/before), `nearest`, `exact`. Duplicate dates
use the last entry for last/exact; nearest ties use the first sorted entry.
Unmatched or overly distant factors are `NaN` with a warning. `gap_days` is signed.

Manifests may be CSV, TSV, TXT/TAB, XLSX, XLS or a pandas DataFrame. Columns are
case insensitive: `filename` (required), `isotope` (unless `half_life`),
`calibration_factor`, `pet_start`, `sub`, `ses`, `half_life`. Relative filenames
resolve against the manifest folder (the current directory for DataFrames).
File manifests preserve string labels such as `01`; use string labels when
constructing a DataFrame too. Each missing factor requires `cal_dates`/`cal_values`
and is recorded in `calibration_factors.tsv`. Invalid rows, including failed
calibration reads, are reported without stopping the remaining rows. A global
`bids_dir` applies only to rows with a subject. Additional options such as
`frame_scheme`, `recording` and `overwrite` are forwarded to processing.

## Tests

```sh
python -m pytest Python/tests
```

Synthetic tests cover arithmetic, onset detection, framing, dates, manifests,
QC and exports. Tests comparing against R are skipped when R or its dependencies
are unavailable locally and run in GitHub CI. Vendor reference data are not
bundled; no independent vendor-level validation of the Python port is claimed.
