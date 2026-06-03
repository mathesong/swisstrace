
# swisstrace

<!-- badges: start -->
<!-- badges: end -->

The goal of swisstrace is to enable automatic calibration and correction of 
Swisstrace `.crv` files, and to facilitate conversion of this data to the BIDS
standard.

## Installation

You can install the development version of swisstrace like so:

``` r
# install.packages("remotes")
remotes::install_github("mathesong/swisstrace")
```

## Usage

`swisstrace_process()` is the main entry point: it takes a raw twilite `.crv`
recording, corrects and calibrates it, and writes the outputs to disk — the
calibrated whole-blood curve for PMOD, a QC plot, and BIDS PET blood files.

The typical call supplies four things: the raw `.crv` file, the per-study calibration
value, the isotope (which sets the decay half-life), and the scan start time.

``` r
library(swisstrace)

swisstrace_process(
  "ABC_P009_D1.crv",           # raw twilite recording
  calibration_factor = 0.425,  # per-study calibration (counts/sec -> kBq/cc)
  isotope = "F18",             # required: not stored in the .crv; sets the half-life
  pet_start = "11:16:30"       # scan start ("time 0") as a clock time, HH:MM:SS
)
```

This writes, under the file's folder (or under a specified `output_folder`):

  - `Corrected_PMOD/ABC_P009_D1_corrected.crv` — the calibrated curve for use with PMOD
  - `Plots/ABC_P009_D1_corrected.png` — a QC plot

Add `sub` (and optionally `ses`) to also emit BIDS PET blood files under
`BIDS/sub-<sub>/[ses-<ses>/]pet/`:

``` r
swisstrace_process("ABC_P009_D1.crv", calibration_factor = 0.425,
                   isotope = "F18", pet_start = "11:16:30",
                   sub = "P009", ses = "01")
```

To write the blood files straight into an existing BIDS dataset, pass `bids_dir`. The
function looks for that subject's (and session's) `_pet.nii.gz` image inside the
dataset and drops the blood `.tsv`/`.json` next to it, named from that image's BIDS
entities so the recording is tied to the scan. If no matching image is found, it falls
back to `bids_dir/sub-<sub>/[ses-<ses>/]pet/`.

``` r
swisstrace_process("ABC_P009_D1.crv", calibration_factor = 0.425,
                   isotope = "F18", pet_start = "11:16:30",
                   sub = "P009", ses = "01", 
                   bids_dir = "/data/my_study")
```

### The scan start time

If `pet_start` is omitted entirely, this time set automatically to 20 seconds 
before the rise to the peak begins. This should produce little to no bias for
radiotracers with long half-lives, but may be non-negligible especially for 15O
studies.

### Output sampling (frames)

By default the corrected curve is returned **as-is**, at the twilite's native
sampling — no resampling is applied. To resample onto a frame scheme instead, pass a
`frame_scheme` data.frame with `width` and `end` columns (both in seconds). Copy the
example below and adjust the widths and ends to taste — here, 1-second frames up to
180 s, then 10-second frames up to 600 s:

``` r
swisstrace_correct(
  "ABC_P009_D1.crv", calibration_factor = 0.425, isotope = "F18",
  frame_scheme = data.frame(width = c(1, 10),
                            end   = c(180, 600))
)
```

This also works with `swisstrace_process()` and `swisstrace_convert_batch()`, which
forward `frame_scheme` (and other options) through to `swisstrace_correct()`.

### Looking up the calibration value

The calibration factor is established periodically for the system. If you keep a
record of those calibrations (their dates and values), `lookup_calibration()` picks the
one that applies to a given measurement — by default the most recent calibration on or before
the study date. You supply the study together with your recorded calibration dates and
values. Identify the study in one of two ways: `filename =` a path to
an existing `.crv` (its acquisition date is read from the recording), or `date =` the
study date directly as `YYYY-MM-DD` or `YYYY/MM/DD` (a `Date`/`POSIXct` or a result from
the functions above also work). Supply exactly one:

``` r
calibrations <- data.frame(
  cal_dates  = c("2025-09-04", "2026-01-15", "2026-02-20"),
  cal_values = c(0.198, 0.41, 0.425)
)

# from a .crv file (date read from the recording)
cal <- lookup_calibration(
  filename = "ABC_P009_D1.crv",
  calibrations$cal_dates,
  calibrations$cal_values
)

# or directly from a study date
cal <- lookup_calibration(
  date = "2026/02/22",
  calibrations$cal_dates,
  calibrations$cal_values
)
cal$calibration_factor    # the factor to pass as `calibration_factor`
```

## Converting many files at once

`swisstrace_convert_batch()` runs the whole pipeline over a list of recordings
described in a manifest — an Excel, CSV, or TSV file (or a data.frame) with one row per
recording. Recognised columns (case-insensitive):

| column | required | meaning |
|---|---|---|
| `filename` | yes | path to the raw `.crv` (absolute, or relative to the manifest) |
| `isotope` | yes\* | isotope name, e.g. `F18` (\*or give `half_life` instead) |
| `pet_start` | no | scan start `"HH:MM:SS"` or seconds; blank = auto-detect |
| `calibration_factor` | no | per-study factor; looked up when absent/blank |
| `sub`, `ses` | no | BIDS labels — enable BIDS output for that row |
| `half_life` | no | half-life in seconds, for isotopes off the built-in table |

``` r
swisstrace_convert_batch("manifest.csv", output_folder = "converted")
```

Each row is handed to `swisstrace_process()`, so the same outputs are written (the
calibrated curve, a QC plot, and BIDS files for rows that carry a `sub`). One bad row
is reported but does not stop the rest; the returned table summarises each row's
status.

If some rows have no `calibration_factor`, pass your calibration record and it is
filled in per study via `lookup_calibration()`. A `calibration_factors.tsv`
(`filename`, `date`, `calibration_factor`, `gap_days`) is then written to the output
folder so you can check what was used:

``` r
swisstrace_convert_batch(
  "manifest.csv",
  output_folder = "converted",
  cal_dates  = c("2025-09-04", "2026-01-15", "2026-02-20"),
  cal_values = c(0.198, 0.41, 0.425),
  bids_dir   = "/data/my_study"     #i.e. drop BIDS blood files into the dataset
)
```

