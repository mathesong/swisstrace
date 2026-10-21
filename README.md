# swisstrace

Correct and calibrate raw Swisstrace **twilite** automatic blood-sampler `.crv`
recordings into whole-blood activity curves, with PMOD, QC plots and BIDS PET blood
exports. Equivalent tools are available in R, Python and MATLAB.

| Folder | Contents | Installation from this checkout |
| --- | --- | --- |
| [R/](R/) | Complete R package: source, metadata, documentation and tests | `remotes::install_local("R")` |
| [Python/](Python/) | Installable Python package, documentation and tests | `python -m pip install ./Python` |
| [MATLAB/](MATLAB/) | MATLAB functions, documentation and tests; no additional toolboxes | `addpath("MATLAB")` |

The repository root is no longer an R package. Existing R installations should
use the `R` subdirectory:

```r
remotes::install_github("CPernet/swisstrace", subdir = "R")
library(swisstrace)
swisstrace_process("recording.crv", calibration_factor = 0.425,
                   isotope = "F18", pet_start = "11:16:30")
```

```python
from swisstrace import swisstrace_process

paths = swisstrace_process("recording.crv", calibration_factor=0.425,
                           isotope="F18", pet_start="11:16:30")
```

```matlab
addpath('MATLAB');
paths = swisstrace_process('recording.crv', 0.425, 'F18', ...
    'pet_start', '11:16:30');
```

See the [R guide](R/README.md), [Python guide](Python/README.md) or
[MATLAB guide](MATLAB/README.md) for batch
conversion, calibration lookup, frame schemes and BIDS output.

All implementations use the coincidence channel (the first counter), subtract
background, correct decay relative to PET start, and apply the calibration factor.
Native sampling is retained unless a frame scheme is supplied. Supply the correct
isotope or an explicit half-life in seconds; it cannot be inferred from the file.
The first corrected sample/frame is zeroed by default. These are whole-blood
curves, without plasma, metabolite or dispersion correction.

## Development and validation

```sh
python -m pip install -e './Python[test,excel]'
python -m pytest Python/tests
Rscript -e 'testthat::test_local("R")'
```

R development commands such as `devtools::load_all()`, `devtools::document()` and
`devtools::install()` run inside `R/` (or receive `"R"` as their package path).
The R package's own source directory is consequently `R/R/`.

Tests use synthetic recordings, including known analytical correction values.
Direct R/Python comparisons run when R plus `tibble` and `jsonlite` are available;
otherwise those tests are skipped. GitHub CI installs both runtimes and checks
numerical agreement, calibration lookup and exported curves. Real recordings and
vendor reference outputs are not included, so the Python port has not been
independently validated against those reference outputs.

MATLAB support is tested on R2024b on Linux; other releases and platforms are
unverified. No R2021a+ compatibility claim is made. Run its tests with:

```matlab
addpath('MATLAB');
results = runtests('MATLAB/tests');
assertSuccess(results);
```

The MATLAB workflow generates reference outputs directly from `R/R/`, then runs
MATLAB tests comparing correction, framing, calibration lookup and exported curves.
For those comparisons locally, first run
`Rscript MATLAB/tests/generate_r_reference.R` from the repository root; otherwise
the two R comparison tests are skipped. Synthetic tests do not replace validation
against real twilite/vendor reference recordings.

Real/vendor validation is a documented follow-up for the MATLAB port, rather than
a blocker to merging the synthetic/R-parity-tested implementation. It is required
before claiming vendor equivalence or relying on the port for production analysis;
see the [validation follow-up](MATLAB/README.md#realvendor-validation-follow-up).
