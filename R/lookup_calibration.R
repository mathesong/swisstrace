# Match a per-study calibration factor to an acquisition by date.
#
# Calibration factors are established periodically (a calibration-vs-date sheet).
# The factor to use for a study is the one from the most recent calibration on or
# before the acquisition date (the manual's "last calibration procedure").

#' Read the acquisition start date from a raw twilite .crv
#'
#' @param file Path to a raw .crv file.
#' @param tz Timezone for interpreting the timestamp (default "UTC").
#' @return A `Date`.
#' @export
swisstrace_date <- function(file, tz = "UTC") {
  stopifnot(file.exists(file))
  ln <- readLines(file, n = 50, warn = FALSE)
  ln <- trimws(ln[nzchar(trimws(ln))])
  p <- NULL
  for (l in ln) {
    v <- suppressWarnings(as.numeric(strsplit(l, "[[:space:]]+")[[1]]))
    if (length(v) >= 6 && all(is.finite(v[1:6]))) { p <- v; break }
  }
  if (is.null(p)) stop("Could not read a timestamp from '", file, "'.")
  as.Date(ISOdatetime(p[1], p[2], p[3], p[4], p[5], p[6], tz = tz))
}

# Read acquisition Date(s) from one or more raw .crv paths.
.dates_from_crv <- function(files, tz = "UTC") {
  missing <- !file.exists(files)
  if (any(missing)) {
    stop("File not found: ",
         paste(sprintf("'%s'", files[missing]), collapse = ", "), ".")
  }
  as.Date(vapply(files, function(f) as.character(swisstrace_date(f, tz)),
                 character(1)))
}

# Coerce a `date` argument (no files) to a vector of acquisition Dates.
.as_study_date <- function(x) {
  if (inherits(x, "Date"))    return(x)
  if (inherits(x, "POSIXct")) return(as.Date(x))
  if (is.list(x) && !is.null(x$date)) return(as.Date(x$date))  # swisstrace_correct() result
  if (is.character(x)) {
    out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
    for (k in seq_along(x)) {
      xk <- trimws(x[k])
      if (grepl("^\\d{4}[-/]\\d{1,2}[-/]\\d{1,2}$", xk)) {
        out[k] <- as.Date(gsub("/", "-", xk))    # YYYY-MM-DD or YYYY/MM/DD literal
      } else {
        stop("Could not interpret '", x[k], "' as a date. ",
             "Use \"YYYY-MM-DD\"/\"YYYY/MM/DD\", or pass a .crv path via `filename`.")
      }
    }
    return(out)
  }
  as.Date(x)
}

#' Look up the calibration factor for a study by date
#'
#' Identify the study either by `filename` (a raw `.crv` whose acquisition date is
#' read from the file) or by `date` (given directly). Supply exactly one of the two.
#'
#' @param cal_dates Vector of calibration dates (Date or coercible).
#' @param cal_values Vector of recorded calibration factors (same length as `cal_dates`).
#' @param filename Path(s) to an existing raw `.crv` file; the acquisition start
#'   date is read from each. Mutually exclusive with `date`. Vectorised.
#' @param date The acquisition date(s), given directly. Accepts a `Date`/`POSIXct`,
#'   a date string ("YYYY-MM-DD" or "YYYY/MM/DD"), or a `swisstrace_correct()` result
#'   list (its `$date` is used). Mutually exclusive with `filename`. Vectorised.
#' @param method How to match: "last" (most recent calibration on or before the study
#'   date; default), "nearest" (closest in either direction), or "exact".
#' @param max_gap Optional maximum allowed gap in days between the study and the
#'   matched calibration; matches beyond it are returned as `NA` with a warning.
#' @param tz Timezone for reading dates from `.crv` files (default "UTC").
#'
#' @return A tibble with one row per study: `study_date`, `calibration_factor`,
#'   `cal_date` (the calibration used), and `gap_days` (study_date - cal_date;
#'   positive = calibration predates the study).
#' @export
lookup_calibration <- function(cal_dates, cal_values, filename = NULL, date = NULL,
                               method = c("last", "nearest", "exact"),
                               max_gap = NULL, tz = "UTC") {
  method <- match.arg(method)
  if (is.null(filename) && is.null(date)) {
    stop("Supply a study via `filename` (a .crv path) or `date`.")
  }
  if (!is.null(filename) && !is.null(date)) {
    stop("Supply only one of `filename` or `date`, not both.")
  }
  study <- if (!is.null(filename)) .dates_from_crv(filename, tz) else .as_study_date(date)
  cal_dates  <- as.Date(cal_dates)
  cal_values <- as.numeric(cal_values)
  if (length(cal_dates) != length(cal_values)) {
    stop("`cal_dates` and `cal_values` must have the same length.")
  }
  if (anyNA(cal_dates)) stop("`cal_dates` contains unparseable dates.")
  ord <- order(cal_dates)
  cal_dates <- cal_dates[ord]; cal_values <- cal_values[ord]

  match_one <- function(d) {
    if (is.na(d)) return(NA_integer_)
    if (method == "exact") {
      hit <- which(cal_dates == d)
      return(if (length(hit)) hit[length(hit)] else NA_integer_)
    }
    if (method == "last") {
      elig <- which(cal_dates <= d)
      return(if (length(elig)) elig[length(elig)] else NA_integer_)
    }
    which.min(abs(as.numeric(cal_dates - d)))  # nearest
  }
  idx <- vapply(study, match_one, integer(1))

  cal_date <- cal_dates[idx]
  gap <- as.numeric(study - cal_date)
  value <- cal_values[idx]
  if (!is.null(max_gap)) {
    too_far <- !is.na(gap) & abs(gap) > max_gap
    value[too_far] <- NA_real_
  }

  if (anyNA(value)) {
    warning(sum(is.na(value)), " of ", length(value),
            " stud", ifelse(length(value) == 1, "y", "ies"),
            " had no calibration factor within range",
            if (method == "last") " (none on or before the study date; try method = 'nearest')" else "",
            ".")
  }

  tibble::tibble(
    study_date         = study,
    calibration_factor = value,
    cal_date           = cal_date,
    gap_days           = gap
  )
}
