# swisstrace for MATLAB

MATLAB equivalents of all seven public [R functions](../R/README.md), ported
from `R/R/`. Supported/tested configuration: **MATLAB R2024b on Linux**,
with no additional toolboxes. Other releases and platforms are unverified;
R2021a+ compatibility is not currently claimed.
R and Python are not needed to use this implementation.

## Installation

Clone or download the repository, then add the MATLAB folder to your path:

```matlab
addpath('/path/to/swisstrace/MATLAB');
% Optional: savepath  % keep this path for future MATLAB sessions
```

From the repository root, `addpath('MATLAB')` is sufficient. The private helper
folder is found automatically; do not add it separately. No pip/uv installation
or compilation is needed. This implementation uses MATLAB tables and datetimes;
GNU Octave compatibility is not claimed.

## Correct and export

```matlab
res = swisstrace_correct('recording.crv', 0.425, 'F18', ...
    'pet_start', '11:16:30');
disp(res.tac);  % time, activity, frame_start, frame_end, frame_dur
fig = swisstrace_qc(res);

paths = swisstrace_process('recording.crv', 0.425, 'F18', ...
    'pet_start', '11:16:30', 'output_folder', 'converted', ...
    'sub', '01', 'ses', '02');
```

`paths` holds the corrected PMOD `.crv`, QC `.png` and BIDS `.tsv`/`.json` paths.
Without `sub`, only PMOD and QC files are written. Default output locations:

- `Corrected_PMOD/<base>_corrected.crv`
- `Plots/<base>_corrected.png`
- `BIDS/sub-01/ses-02/pet/sub-01_ses-02_recording-autosampler_blood.tsv` and `.json`

Pass `bids_dir` to write into an existing dataset. Matching PET image filenames
supply the output directory and entities; multiple matches warn and select the
first sorted path, following R. If no image matches, the subject/session PET
folder is used. `ses` or `bids_dir` requires a subject. Keep labels as strings to
preserve leading zeros. `recording` defaults to `'autosampler'`.
`overwrite=false` preserves each existing output individually.

## Functions and arguments

| Function | MATLAB return |
| --- | --- |
| `swisstrace_correct(file, factor, isotope, ...)` | Struct; `tac` and `raw` are tables |
| `swisstrace_process(file, factor, isotope, ...)` | Struct of output paths |
| `swisstrace_qc(res, ...)` | Figure handle; use `close(fig)` when finished |
| `swisstrace_date(file, tz)` | Timezone-free datetime at midnight |
| `lookup_calibration(dates, values, ...)` | Table of factors, matched dates and signed gaps |
| `swisstrace_convert_batch(manifest, ...)` | Table with per-row status and error message |
| `assert_raw_crv(file)` | `true`, or an error for non-raw files |

Required correction inputs are positional; optional inputs use name/value pairs.
Use `[]` for R's `NULL`. Use `res.tac`/`res.background` to access result fields.
Function help is available through `help swisstrace_correct`, etc.

Correction options retain the R defaults:

| Option | Default | Meaning |
| --- | --- | --- |
| `pet_start` | `[]` | Auto-detect time zero, or supply seconds/datetime/clock time |
| `half_life` | `[]` | Override isotope lookup, in seconds |
| `lead` | `20` | Seconds before detected onset to place time zero |
| `frame_scheme` | `[]` | Native samples; or a width/end table or `[width end]` matrix |
| `zero_first_frame` | `true` | Force first sample/frame to zero |
| `baseline_k` | `3` | Detection threshold in scaled median absolute deviations |
| `min_run` | `5` | Consecutive elevated samples required |
| `baseline_init` | `10` | Initial baseline duration, seconds |
| `baseline_min` | `20` | Warn when detected onset has less baseline, seconds |
| `tz` | `'UTC'` | Timezone for recorded wall-clock timestamps |

Isotope names are case/punctuation insensitive: F18, C11, N13, O15, Ga68, Cu62,
Zr89, I124 and Rb82. To use a direct half-life, pass an empty isotope:

```matlab
res = swisstrace_correct('recording.crv', .425, [], 'half_life', 6586.2);
res = swisstrace_correct('recording.crv', .425, 'F18', ...
    'pet_start', 40, 'frame_scheme', [1 180; 10 600], ...
    'zero_first_frame', false);
```

The first counter is the coincidence channel; singles are retained only in
`res.raw`. Background subtraction, decay correction and calibration follow R.
Automatic detection uses R's `1.4826 * median(abs(x - median(x)))`, falling back
to sample standard deviation, then 1. Explicit PET start uses all preceding
samples for background; automatic mode uses the accumulated clean baseline.

Clock strings (`HH:MM[:SS[.fraction]]`) use the acquisition's local date. For a
scan after midnight, pass elapsed seconds or a full datetime. A datetime without
a timezone is interpreted in `tz`. As in R, acquisition `date` and datetime
calibration inputs use the UTC calendar date; acquisition datetimes retain `tz`.

Frames average samples in left-closed, right-open intervals. Empty frames are
`NaN`; exported missing values are `NA` for PMOD and `n/a` for BIDS. Native
sample times and fractional seconds are retained. Time metadata and all result
fields use the R names, including `lambda`, `t0_seconds` and `n_background`.

## Calibration and batch conversion

```matlab
cal_dates = ["2025-09-04", "2026-01-15", "2026-02-20"];
cal_values = [.198, .41, .425];
cal = lookup_calibration(cal_dates, cal_values, ...
    'filename', 'recording.crv', 'method', 'last', 'max_gap', 90);
% Alternatively: 'date', '2026-02-22' (do not supply both date and filename).

status = swisstrace_convert_batch('manifest.xlsx', ...
    'output_folder', 'converted', 'cal_dates', cal_dates, 'cal_values', cal_values);
disp(status(:, {'filename','status','message'}));
```

Lookup supports vector inputs, date strings (`YYYY-MM-DD` or `YYYY/MM/DD`),
datetime values, and a correction result as `date`. Methods are `last`, `nearest`
and `exact`. Last/exact choose the last duplicate date; nearest ties choose the
first sorted entry. An unmatched/overly distant factor becomes `NaN` with a
warning; matched date and signed `gap_days` are retained for max-gap exclusions.

Manifests may be MATLAB tables or CSV/TSV/TXT/TAB/XLSX/XLS files. Recognised
columns (case insensitive): `filename`, `isotope`, `half_life`, `pet_start`,
`calibration_factor`, `sub`, `ses`. Filenames are required; isotope or half-life
is required per row. Files are imported as text to preserve labels like `01`.
Relative recording paths resolve from the manifest folder (current directory
for table input). Missing calibration factors use the supplied log and are
recorded in `calibration_factors.tsv`. Additional correction/export options
apply to all rows; manifest columns override the corresponding shared options.
A global `bids_dir` is used only for rows with a subject.

## Differences and validation

The MATLAB port rejects nonpositive/nonfinite factors or half-lives, invalid
frame schemes/dates, out-of-order timestamps, and starts after recording end.
It handles short flat recordings and isolates per-row calibration lookup errors
in batches. QC returns a figure rather than R's invisible result. BIDS time
descriptions distinguish native sample times from frame midpoints. These are
intentional interface/validation differences; valid-input correction follows R.

```matlab
addpath('MATLAB');
results = runtests('MATLAB/tests');
assertSuccess(results);
```

Tests cover synthetic analytical signals, all isotopes, onset detection, frame
boundaries/missing frames, dates, calibration matching, plots, exports and batch
formats. To enable the two direct R comparison tests locally, run this first
from the repository root (requires R with `tibble` and `jsonlite`):

```sh
Rscript MATLAB/tests/generate_r_reference.R
```

The MATLAB GitHub workflow generates those references from the actual R code and
compares native/framed/automatic/C11 outputs and exports. Without generated
references, only those two comparison tests are skipped. Generated fixtures are
ignored by Git. Real/vendor reference recordings are not bundled, so independent
validation against those data remains necessary for the MATLAB port.

### Real/vendor validation follow-up

The analytical tests and direct R comparisons establish implementation parity on
synthetic inputs. This is sufficient to merge the MATLAB port with the limitation
above; it does not establish independent agreement with vendor outputs. Real/vendor
validation remains a follow-up and is required before claiming vendor equivalence
or relying on this port for production analysis.

To complete that follow-up:

- Use representative raw twilite recordings with matched vendor/reference corrected
  curves, known calibration factors, isotope/half-life and PET start. Keep restricted
  recordings outside Git; commit shareable fixtures only with permission.
- Compare MATLAB with both R and the independent reference, checking activity,
  time alignment, units, native/framed sampling and automatic versus supplied starts.
- Record MATLAB release/platform, input provenance, settings and maximum absolute
  activity error. Define the acceptable tolerance for each reference before comparing;
  the historical R error below 0.06 kBq/cc is not a validated MATLAB tolerance.
- Investigate discrepancies, add shareable regression coverage where possible and
  document the results and remaining limitations before extending validation claims.
