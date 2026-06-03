# Batch-convert many twilite .crv files listed in a manifest spreadsheet.
#
# The manifest (an .xlsx/.xls/.csv/.tsv file, or a data.frame) has one row per
# recording. Recognised columns (case-insensitive):
#   filename            path to the raw .crv (absolute, or relative to the manifest)
#   isotope             isotope name, e.g. "F18"  (or give half_life instead)
#   pet_start           scan start "HH:MM:SS" or seconds; blank = auto-detect (optional)
#   calibration_factor  per-study factor (optional; looked up when absent/blank)
#   sub, ses            BIDS labels (optional, enable BIDS output for that row)
#   half_life           half-life in seconds (optional; for isotopes off the table)
#
# Each row is handed to swisstrace_process(). When a row has no calibration_factor,
# it is filled with lookup_calibration() from the `cal_dates`/`cal_values` log, and a
# `calibration_factors.tsv` (filename, date, calibration_factor, gap_days) is written
# to the output folder. One bad row does not abort the batch.

# Read a manifest (path or data.frame) into a data.frame with lower-cased column names.
.read_manifest <- function(x) {
  if (is.data.frame(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    if (!is.character(x) || length(x) != 1L) {
      stop("`manifest` must be a single file path or a data.frame.")
    }
    if (!file.exists(x)) stop("Manifest file not found: '", x, "'.")
    ext <- tolower(tools::file_ext(x))
    df <- switch(ext,
      xlsx = ,
      xls = {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("Reading an Excel manifest needs the 'readxl' package. Install it with ",
               "install.packages(\"readxl\"), or save the manifest as .csv or .tsv.")
        }
        as.data.frame(readxl::read_excel(x, col_types = "text"), stringsAsFactors = FALSE)
      },
      tsv = ,
      tab = ,
      txt = utils::read.delim(x, colClasses = "character",
                              na.strings = c("NA", ""), check.names = TRUE),
      csv = utils::read.csv(x, colClasses = "character",
                            na.strings = c("NA", ""), check.names = TRUE),
      stop("Unsupported manifest type '.", ext, "'. Use .xlsx, .xls, .csv, or .tsv.")
    )
  }
  names(df) <- tolower(trimws(names(df)))
  df
}

# Turn a manifest pet_start cell into a value for swisstrace_correct(): NULL (auto),
# a number of seconds, or a clock-time string left as-is for "HH:MM:SS" parsing.
.parse_pet_start <- function(x) {
  if (length(x) != 1L || is.na(x)) return(NULL)
  if (is.character(x)) {
    x <- trimws(x)
    if (!nzchar(x)) return(NULL)
    if (grepl("^[0-9]*\\.?[0-9]+$", x)) return(as.numeric(x))  # seconds given as a number
    return(x)                                                  # "HH:MM:SS"
  }
  x                                                            # already numeric
}

#' Batch-convert twilite .crv files listed in a manifest
#'
#' Reads a manifest spreadsheet describing many recordings and runs
#' [swisstrace_process()] on each. See the columns below.
#'
#' @param manifest Path to an `.xlsx`/`.xls`/`.csv`/`.tsv` file, or a data.frame, with
#'   one row per recording. Recognised columns (case-insensitive): `filename`
#'   (required), `isotope` (required unless `half_life` is given), `pet_start`
#'   (optional - `"HH:MM:SS"` or seconds; blank auto-detects), `calibration_factor`
#'   (optional - looked up when absent/blank), `sub`/`ses` (optional BIDS labels), and
#'   `half_life` (optional, seconds). Relative `filename`s are resolved against the
#'   manifest's folder.
#' @param output_folder Folder for the outputs (and the calibration log). Defaults to
#'   the manifest's folder, or the working directory if `manifest` is a data.frame.
#' @param cal_dates,cal_values Your record of calibrations (dates and factors), used to
#'   fill in `calibration_factor` for any row that lacks one, via [lookup_calibration()].
#'   Required only if some rows need a lookup.
#' @param bids_dir Optional BIDS dataset root passed through to [swisstrace_process()]
#'   (used for rows that carry a `sub`).
#' @param cal_method,cal_max_gap Passed to [lookup_calibration()] (`method`, `max_gap`).
#' @param ... Further arguments forwarded to [swisstrace_process()] for every row
#'   (e.g. `recording`, `overwrite`, `lead`, `frame_scheme`).
#'
#' @return (invisibly) a tibble with one row per manifest entry: `filename`, `isotope`,
#'   `calibration_factor`, `status` ("ok"/"error"), and `message` (the error, if any).
#' @export
swisstrace_convert_batch <- function(manifest,
                                     output_folder = NULL,
                                     cal_dates = NULL,
                                     cal_values = NULL,
                                     bids_dir = NULL,
                                     cal_method = c("last", "nearest", "exact"),
                                     cal_max_gap = NULL,
                                     ...) {
  cal_method <- match.arg(cal_method)
  man <- .read_manifest(manifest)

  miss <- setdiff("filename", names(man))
  if (length(miss)) {
    stop("Manifest is missing required column(s): ", paste(miss, collapse = ", "),
         ". Found: ", paste(names(man), collapse = ", "), ".")
  }
  n <- nrow(man)

  manifest_dir <- if (is.character(manifest) && length(manifest) == 1L &&
                      !is.data.frame(manifest) && file.exists(manifest)) {
    dirname(normalizePath(manifest))
  } else {
    getwd()
  }
  output_folder <- output_folder %||% manifest_dir
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

  # resolve each filename: as given, else relative to the manifest's folder
  resolve <- function(f) {
    f <- as.character(f)
    if (is.na(f) || !nzchar(f)) return(NA_character_)
    if (file.exists(f)) return(f)
    alt <- file.path(manifest_dir, f)
    if (file.exists(alt)) alt else f
  }
  paths <- vapply(man$filename, resolve, character(1))

  ## --- resolve calibration factors -----------------------------------------
  cf <- if ("calibration_factor" %in% names(man)) {
    suppressWarnings(as.numeric(man$calibration_factor))
  } else {
    rep(NA_real_, n)
  }
  need <- is.na(cf)
  cal_log_path <- NULL
  if (any(need)) {
    if (is.null(cal_dates) || is.null(cal_values)) {
      stop(sum(need), " of ", n, " row(s) have no `calibration_factor`, and no ",
           "`cal_dates`/`cal_values` were supplied to look it up. Either give the ",
           "factor in the manifest, or pass your calibration log via `cal_dates` and ",
           "`cal_values`.")
    }
    look <- lookup_calibration(cal_dates, cal_values, filename = paths[need],
                               method = cal_method, max_gap = cal_max_gap)
    cf[need] <- look$calibration_factor

    cal_log <- tibble::tibble(
      filename           = man$filename[need],
      date               = look$study_date,
      calibration_factor = look$calibration_factor,
      gap_days           = look$gap_days
    )
    cal_log_path <- file.path(output_folder, "calibration_factors.tsv")
    utils::write.table(cal_log, cal_log_path, sep = "\t",
                       row.names = FALSE, quote = FALSE)
  }

  ## --- process each row -----------------------------------------------------
  cell <- function(col, i) if (col %in% names(man)) man[[col]][i] else NA
  label <- function(v) {       # a BIDS label, or NULL when blank/missing
    v <- as.character(v)
    if (length(v) != 1L || is.na(v) || !nzchar(trimws(v))) NULL else trimws(v)
  }
  dots <- list(...)
  status  <- character(n)
  message_ <- rep(NA_character_, n)

  for (i in seq_len(n)) {
    out <- tryCatch({
      if (is.na(cf[i])) {
        stop("no calibration factor (absent from the manifest and no match found in ",
             "the calibration log).")
      }
      hl <- suppressWarnings(as.numeric(cell("half_life", i)))
      args <- list(
        filename           = paths[i],
        calibration_factor = cf[i],
        isotope            = as.character(cell("isotope", i)),
        pet_start          = .parse_pet_start(cell("pet_start", i)),
        output_folder      = output_folder,
        sub                = label(cell("sub", i)),
        ses                = label(cell("ses", i)),
        bids_dir           = bids_dir
      )
      if (!is.na(hl)) args$half_life <- hl
      suppressMessages(do.call(swisstrace_process, c(args, dots)))
      "ok"
    }, error = function(e) structure("error", msg = conditionMessage(e)))
    status[i] <- out
    if (identical(out, "ok")) {
      message(sprintf("[%d/%d] %s: ok", i, n, basename(paths[i])))
    } else {
      message_[i] <- attr(out, "msg")
      message(sprintf("[%d/%d] %s: ERROR - %s", i, n, basename(paths[i]), message_[i]))
    }
  }

  n_ok <- sum(status == "ok")
  message(sprintf("\nConverted %d/%d file(s).%s", n_ok, n,
                  if (!is.null(cal_log_path))
                    paste0(" Calibration log: ", cal_log_path) else ""))

  invisible(tibble::tibble(
    filename           = man$filename,
    isotope            = if ("isotope" %in% names(man)) as.character(man$isotope) else NA_character_,
    calibration_factor = cf,
    status             = status,
    message            = message_
  ))
}
