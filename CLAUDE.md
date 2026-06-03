# swisstrace

An R package that corrects and calibrates raw Swisstrace **twilite**
automatic-blood-sampler `.crv` files (whole-blood radioactivity over time) into
calibrated whole-blood activity input functions (AIFs, kBq/cc). The goal is a generic,
automated, scriptable correction so this step can be done with ordinary tooling in R.

## Layout

This repo is the **R package itself** (standard layout):

- `R/` — package source.
  - `swisstrace_correct.R` — the core correction (also defines `assert_raw_crv()`).
  - `lookup_calibration.R` — `lookup_calibration()` + `swisstrace_date()`.
  - `swisstrace_qc.R` — `swisstrace_qc()`.
  - `swisstrace_process.R` — `swisstrace_process()`.
  - `swisstrace_convert_batch.R` — `swisstrace_convert_batch()` (manifest-driven batch).
- `man/` — roxygen-generated `.Rd` docs (do not hand-edit; run `devtools::document()`).
- `DESCRIPTION`, `NAMESPACE` — package metadata (R 4.x; Imports: graphics,
  grDevices, jsonlite, stats, tibble, tools, utils; Suggests: readxl for Excel
  manifests, testthat).

Build/iterate with `devtools::load_all()` / `devtools::document()` /
`devtools::install()`. Real `.crv` data and reference outputs live outside this repo,
in a separate data/analysis workspace.

## CRV format

Whitespace-delimited, no header. Each row:
```
YYYY M D H M S.sss   <Counter_1>   <Counter_2>   <Counter_3>
```
- **Counter_1 = coincidences (counts/sec) — the signal used for correction.**
- Counter_2, Counter_3 = the two singles channels (LYSO crystals); **not used**.

`assert_raw_crv()` rejects already-corrected/processed `.crv` files up front (corrected
files begin with a `Corrected_&_calibrated ... value` header). It is called by both
`swisstrace_correct()` and `swisstrace_process()`.

## The correction (reverse-engineered, validated to <0.06 kBq/cc against reference outputs)

```
corrected(t) = (coincidence(t) − background) × exp(λ·(t − t0)) × calibration_factor
```
1. **t0 (PET scan start / "time 0")** is *not* stored in the `.crv`, but is usually
   recorded for the scan and worth supplying via `pet_start` (it sets the decay
   reference — small effect for F-18, larger for short half-lives). `pet_start` accepts
   a clock time `"HH:MM:SS"` (taken on the recording's date), a POSIXct, or seconds
   after the recording start. If NULL it is estimated from the curve by forward
   accumulation of the bolus rise onset (`lead`, `baseline_k`, `min_run`,
   `baseline_init`, `baseline_min` control the detection).
2. **Background** = mean coincidence over the baseline before t0. Subtracted from
   every sample.
3. **Decay**: `λ = ln2 / half_life`; F-18 = 6586.2 s. Referenced to t0. Pass
   `isotope` (e.g. "F18", "C11", "O15", "Ga68") or a raw `half_life`. Branching-ratio
   correction is left OFF (folded into the calibration factor).
4. **Calibration factor**: per-study multiplier, from a calibration-vs-date sheet.
5. **Cut** before t0, then **resample** (`frame_scheme`: default 1 s frames to 180 s,
   then 10 s to 600 s; first frame forced to 0 when `zero_first_frame = TRUE`).
   Output = mid-time (s) + activity (kBq/cc).

`isotope` is **required** (a wrong one silently corrupts the decay correction) unless
`half_life` is given directly.

## Usage

```r
res <- swisstrace::swisstrace_correct("ASPC0243_P009_D1.crv",
                                       calibration_factor = 0.425, isotope = "F18")
res$tac   # tidy tibble: time (mid, s), activity (kBq/cc), frame_start/end/dur
# alongside: $date, $pet_start, $background, $half_life, $lambda, $raw, ...
```

Validate against a reference corrected output (should match to ~0.05 kBq/cc):
```r
res <- swisstrace::swisstrace_correct("ASPC0243_P009_D1.crv",
         calibration_factor = 0.425, isotope = "F18", pet_start = 54.5)
ref <- read.table("ASPC0243_P009_260210_corrected.crv", skip = 1, sep = "\t")$V2
max(abs(res$tac$activity - ref))
```

## Exported functions

1. **`swisstrace_correct(file, calibration_factor, isotope, ...)`** — returns a list
   whose `$tac` is the tidy corrected curve, with acquisition date/time + correction
   metadata as sibling elements.
2. **`lookup_calibration(cal_dates, cal_values, filename, date, method, max_gap, ...)`**
   — match a calibration factor to a study by date. Identify the study with exactly one
   of `filename` (a `.crv` path, date read from the file) or `date` (a `Date`/`POSIXct`,
   a "YYYY-MM-DD"/"YYYY/MM/DD" string, or a `swisstrace_correct()` result). Supplying
   neither or both errors. `method`: "last" (default), "nearest", "exact". Vectorised;
   returns a tibble (`study_date`, `calibration_factor`, `cal_date`, `gap_days`).
3. **`swisstrace_date(file, tz)`** — read the acquisition start date from a raw `.crv`.
4. **`swisstrace_process(filename, calibration_factor, isotope, ..., output_folder,
   sub, ses, bids_dir)`** — corrects `filename` (forwarding `...` to
   `swisstrace_correct`) and writes a corrected `.crv` curve, a `<base>_corrected.png`
   QC plot (`Plots/`), and — when `sub` is given — BIDS PET blood `.tsv`/`.json` under
   `BIDS/sub-XX[/ses-XX]/pet/`. `ses` without `sub` errors. If `bids_dir` (a BIDS root)
   is given, the blood files go there instead: beside the subject's
   `sub-XX[_ses-XX]..._pet.nii[.gz]` image when found (named from that image's
   entities), else into `bids_dir/sub-XX[/ses-XX]/pet/`; `bids_dir` requires `sub`.
   BIDS json: PlasmaAvail/MetaboliteAvail false, WholeBloodAvail true,
   DispersionCorrected false.
5. **`swisstrace_qc(res, ...)`** — QC plot of a correction result (raw vs corrected
   traces); returns `res` invisibly.
6. **`swisstrace_convert_batch(manifest, output_folder, cal_dates, cal_values,
   bids_dir, ...)`** — runs `swisstrace_process()` over every row of a manifest
   (`.xlsx`/`.xls`/`.csv`/`.tsv` or a data.frame). Columns (case-insensitive):
   `filename` (req), `isotope` (req unless `half_life`), `pet_start`,
   `calibration_factor`, `sub`/`ses`, `half_life`. Missing factors are filled via
   `lookup_calibration(..., filename=)` from `cal_dates`/`cal_values`, and a
   `calibration_factors.tsv` (filename, date, calibration_factor, gap_days) is written.
   Per-row errors are caught (one bad row doesn't abort); returns a status tibble.

## Conventions / notes

- R, tidyverse-friendly (returns tibbles). Call functions with `swisstrace::`.
- After editing any `R/*.R` roxygen, run `devtools::document()` to regenerate
  `man/` and `NAMESPACE` — never hand-edit those.
- The corrected *values* are nearly independent of t0 (background sits on the flat
  baseline; a small decay-reference shift scales everything <0.2 %). t0's real role is
  the AIF time origin, normally re-aligned via a delay term in kinetic modelling.
- Calibration sheets are institution-specific; their reader is not part of this
  package. `lookup_calibration()` takes the parsed dates/values as plain vectors.

## Next

Batch conversion over many measurements is implemented (`swisstrace_convert_batch()`).
