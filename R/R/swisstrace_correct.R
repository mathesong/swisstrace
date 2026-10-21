# Correction & calibration of Swisstrace twilite blood-sampling .crv files.
#
#   1. The signal is column 1 of the .crv (coincidence counts/sec). Columns 2 & 3
#      are the two singles channels and are not used in the correction.
#   2. A reference time "time 0" (the PET scan start) splits the recording: all
#      samples before it are averaged to give the background, which is subtracted.
#   3. Decay correction exp(lambda * (t - t0)) is applied relative to time 0.
#   4. The result is multiplied by the per-study calibration factor -> kBq/cc.
#   5. The corrected curve is cut before time 0 and resampled onto a frame scheme
#      (default: 1 s frames to 180 s, then 10 s frames to 600 s).
#
# If you know the real PET start time, pass it via `pet_start` (exact). Otherwise it
# is estimated by forward accumulation: seed the baseline from the first
# `baseline_init` s, then walk forward adding each sample that is consistent with the
# baseline (updating its median + MAD) until a run of `min_run` samples rises above
# (baseline median + `baseline_k` * MAD). That run's start is the bolus onset, and
# the accumulated samples are the background. Time 0 is then placed `lead` s before
# the onset (default 20, per swisstrace's guidance to start the twilite >=20 s before
# injection; reduced automatically if less data is available). Defining it by where
# the curve starts to rise makes it robust to injection speed (a slow infusion simply
# rises more gradually). The corrected *values* are almost independent of time 0 (the
# background sits on flat baseline and a small decay-reference shift scales everything
# by <0.2 %); its main role is to set the AIF time origin, which kinetic modelling
# normally re-aligns via a delay term.

# Physical half-lives in seconds.
.tracer_half_lives <- c(
  F18  = 6586.2,    # 109.77 min
  C11  = 1223.4,    # 20.39 min
  N13  = 597.9,     # 9.965 min
  O15  = 122.24,    # 2.037 min
  Ga68 = 4057.74,   # 67.63 min
  Cu62 = 584.4,     # 9.74 min
  Zr89 = 282240,    # 78.4 h
  I124 = 360806.4,  # 4.176 d
  Rb82 = 76.4       # 1.273 min
)

#' Assert that a file is a raw twilite .crv (not a corrected/processed one)
#'
#' Raw recordings begin with a timestamp + counts row; corrected files begin with a
#' `Corrected_&_calibrated ... value` header line. Stops with an informative message
#' if the file looks corrected/non-raw.
#' @param file Path to check.
#' @export
assert_raw_crv <- function(file) {
  stopifnot(file.exists(file))
  ln <- readLines(file, n = 5, warn = FALSE)
  ln <- trimws(ln[nzchar(trimws(ln))])
  if (!length(ln)) stop("File '", file, "' is empty.")
  first <- ln[1]
  looks_corrected <- grepl("corrected|kBq/cc|value\\[|time\\[seconds\\]",
                           first, ignore.case = TRUE)
  numeric_start <- suppressWarnings(!is.na(as.numeric(strsplit(first, "[[:space:]]+")[[1]][1]))) # nolint
  if (looks_corrected || !numeric_start) {
    stop("'", basename(file), "' looks like a corrected/processed .crv, not a raw ",
         "twilite recording.\n",
         "  Expected raw rows of 'YYYY M D H M S  coincidences  singles1  singles2'.\n",
         "  First line was: ", substr(first, 1, 70))
  }
  invisible(TRUE)
}

#' Correct and calibrate a twilite .crv file
#'
#' @param file Path to a raw twilite .crv file.
#' @param calibration_factor Per-study calibration multiplier (counts/sec -> kBq/cc).
#' @param isotope Isotope name for the half-life lookup. REQUIRED: it is not stored
#'   in the `.crv`, and a wrong isotope silently corrupts the decay correction.
#'   E.g. "F18", "C11", "O15", "Ga68". Needed unless `half_life` is supplied directly.
#' @param pet_start PET scan start ("time 0") - the last of the main inputs, alongside
#'   `file`, `calibration_factor` and `isotope`. Usually worth supplying: it sets the
#'   reference for the decay correction (the effect is small for long-lived isotopes
#'   such as F-18 but grows for shorter half-lives). Accepts a clock time as a string
#'   `"HH:MM:SS"` (e.g. `"11:16:30"`, taken on the recording's own date) - the typical
#'   way the scan start is noted - or a POSIXct, or a number of seconds after the start
#'   of the twilite recording. If NULL (the default), time 0 is estimated from the
#'   curve (see `lead`).
#' @param half_life Optional half-life in seconds, used in place of the `isotope`
#'   lookup (for isotopes not in the built-in table).
#' @param lead Seconds to place time 0 before the detected rise onset (only used
#'   when `pet_start` is NULL). Default 20, following swisstrace's guidance to start
#'   the twilite >=20 s before injection. Automatically reduced if fewer than `lead`
#'   seconds were recorded before the onset (it keeps whatever is available). The
#'   lead actually applied is returned as `$lead`.
#' @param frame_scheme Optional data.frame with columns `width` and `end` (seconds)
#'   defining the output resampling. If `NULL` (the default) the corrected samples are
#'   returned as-is at their native sampling, without any resampling. To resample,
#'   pass e.g. `data.frame(width = c(1, 10), end = c(180, 600))` (1 s frames to 180 s,
#'   then 10 s frames to 600 s) and modify to taste.
#' @param zero_first_frame If TRUE (default), the first output frame
#'   is set to exactly 0.
#' @param baseline_k Detection threshold, in robust SDs above the background mean. Default is 3.
#' @param min_run Number of consecutive samples that must exceed the threshold to
#'   mark the rise onset (guards against isolated noise spikes). Default 5.
#' @param baseline_init Seconds at the start of the recording used to seed the
#'   baseline median + MAD before forward accumulation (auto-detect mode). Default 10.
#' @param baseline_min Warn if the rise is detected with less than this many seconds
#'   of accumulated baseline (i.e. the twilite was started too late). Default 20.
#' @param tz Timezone used to interpret the recorded timestamps (default "UTC";
#'   only affects the printed clock time, not any computed quantity).
#' @param ... Reserved for future options; currently ignored.
#'
#' @return A list with:
#'   \item{tac}{tibble: `time` (s; sample time, or frame mid-time when resampled),
#'     `activity` (kBq/cc), and `frame_start`, `frame_end`, `frame_dur`.}
#'   \item{calibration_factor, background, isotope, half_life, lambda}{correction
#'     parameters used.}
#'   \item{date, start_time, acq_start, acq_end, acq_duration}{acquisition timing.}
#'   \item{pet_start, t0_seconds, t0_detected, lead}{the time-0 used and how.}
#'   \item{n_raw, n_background, file, frame_scheme, raw}{provenance, plus the raw
#'     samples as a tibble (`time`, `coincidence`, `singles1`, `singles2`).}
#'
#' @examples
#' \dontrun{
#' # Default: corrected samples returned as-is, at their native sampling.
#' res <- swisstrace_correct("ASPC0243_P009_D1.crv",
#'                           calibration_factor = 0.425, isotope = "F18")
#' res$tac
#'
#' # Resample onto frames: 1 s frames to 180 s, then 10 s frames to 600 s.
#' # Copy this data.frame and modify the widths/ends to taste.
#' res <- swisstrace_correct("ASPC0243_P009_D1.crv",
#'                           calibration_factor = 0.425, isotope = "F18",
#'                           frame_scheme = data.frame(width = c(1, 10),
#'                                                     end   = c(180, 600)))
#' res$tac
#' }
#' @export
swisstrace_correct <- function(file,
                                calibration_factor,
                                isotope,
                                pet_start = NULL,
                                half_life = NULL,
                                lead = 20,
                                frame_scheme = NULL,
                                zero_first_frame = TRUE,
                                baseline_k = 3,
                                min_run = 5,
                                baseline_init = 10,
                                baseline_min = 20,
                                tz = "UTC",
                                ...) {

  stopifnot(file.exists(file))
  assert_raw_crv(file)   # reject corrected/processed .crv files up front
  if (missing(calibration_factor) || !is.numeric(calibration_factor) ||
      length(calibration_factor) != 1) {
    stop("`calibration_factor` must be a single number.")
  }

  ## --- half-life / decay constant ------------------------------------------
  if (is.null(half_life)) {
    if (missing(isotope) || length(isotope) != 1 || is.na(isotope) || !nzchar(isotope)) {
      stop("`isotope` is required: it is not recorded in the .crv, so it must be ",
           "given (e.g. isotope = \"F18\") - a wrong isotope silently corrupts the ",
           "decay correction. Alternatively supply `half_life` (seconds) directly. ",
           "Known isotopes: ", paste(names(.tracer_half_lives), collapse = ", "), ".")
    }
    key <- toupper(gsub("[^A-Za-z0-9]", "", isotope))
    names_lut <- toupper(names(.tracer_half_lives))
    if (!key %in% names_lut) {
      stop("Unknown isotope '", isotope, "'. Supply `half_life` (seconds), or use one of: ",
           paste(names(.tracer_half_lives), collapse = ", "))
    }
    half_life <- unname(.tracer_half_lives[match(key, names_lut)])
  } else if (missing(isotope)) {
    isotope <- NA_character_
  }
  lambda <- log(2) / half_life

  ## --- read raw data --------------------------------------------------------
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(trimws(lines), "[[:space:]]+")
  parts <- parts[vapply(parts, length, integer(1)) >= 9L]
  if (!length(parts)) stop("No parseable data rows found in '", file, "'.")
  m <- vapply(parts, function(p) as.numeric(p[1:9]), numeric(9))
  m <- t(m)

  ts <- ISOdatetime(m[, 1], m[, 2], m[, 3], m[, 4], m[, 5], m[, 6], tz = tz)
  secs <- as.numeric(difftime(ts, ts[1], units = "secs"))
  coinc    <- m[, 7]
  singles1 <- m[, 8]
  singles2 <- m[, 9]

  ## --- determine time 0 -----------------------------------------------------
  robust_sd <- function(x) {
    s <- stats::mad(x)
    if (!is.finite(s) || s == 0) s <- stats::sd(x)
    if (!is.finite(s) || s == 0) s <- 1
    s
  }
  onset_i <- NA_integer_
  bg_idx  <- NULL   # indices of the samples forming the background (auto mode)
  n <- length(coinc)
  if (!is.null(pet_start)) {
    if (inherits(pet_start, "POSIXct")) {
      t0_sec <- as.numeric(difftime(pet_start, ts[1], units = "secs"))
    } else if (is.character(pet_start)) {
      # A clock time "HH:MM:SS" (or "HH:MM") on the acquisition date - the usual way
      # the scan start is recorded. Combined with the recording's own date. Validate
      # the shape first: as.POSIXct() would otherwise silently accept a trailing-junk
      # string as midnight.
      clk <- trimws(pet_start)
      if (!grepl("^[0-9]{1,2}:[0-9]{2}(:[0-9]{2}([.][0-9]+)?)?$", clk)) {
        stop("`pet_start` = \"", pet_start, "\" is not a valid clock time. Use ",
             "\"HH:MM:SS\" (e.g. \"11:16:30\"), a number of seconds after the ",
             "recording start, or a POSIXct.")
      }
      clock <- as.POSIXct(paste(as.Date(ts[1], tz = tz), clk), tz = tz)
      t0_sec <- as.numeric(difftime(clock, ts[1], units = "secs"))
    } else {
      t0_sec <- as.numeric(pet_start)
    }
    t0_detected <- FALSE
    lead_used <- NA_real_
  } else {
    # Forward accumulation: seed the baseline from the first `baseline_init` s, then
    # walk forward. Each sample that is consistent with the baseline is added to it
    # (median + MAD updated); the first sample beginning a run of `min_run` samples
    # above (baseline median + baseline_k * MAD) is the rise onset = time 0.
    base_idx <- which(secs < baseline_init)
    if (length(base_idx) < 3L) base_idx <- seq_len(min(3L, n))
    for (i in seq.int(max(base_idx) + 1L, n)) {
      thr <- stats::median(coinc[base_idx]) + baseline_k * robust_sd(coinc[base_idx])
      j <- min(i + min_run - 1L, n)
      if (coinc[i] > thr && (j - i + 1L) >= min_run && all(coinc[i:j] > thr)) {
        onset_i <- i
        break
      }
      if (coinc[i] <= thr) base_idx <- c(base_idx, i)  # accumulate clean samples only
    }
    if (is.na(onset_i)) {
      warning("No sustained rise above background detected; ",
              "using start of recording as time 0.")
      t0_sec <- 0
      lead_used <- NA_real_
    } else {
      onset_sec <- secs[onset_i]
      if (onset_sec < baseline_min) {
        warning("Rise detected after only ", round(onset_sec, 1),
                " s of baseline (< baseline_min = ", baseline_min,
                " s); the background estimate may be unreliable. ",
                "Start the twilite earlier, or supply `pet_start`.")
      }
      # place time 0 `lead` s before the onset, but never before the recording
      lead_used <- min(lead, onset_sec)
      if (lead_used < lead) {
        message("Only ", round(onset_sec, 1), " s recorded before the rise; ",
                "using a ", round(lead_used, 1), " s lead instead of ", lead, " s.")
      }
      t0_sec <- onset_sec - lead_used
    }
    bg_idx <- base_idx
    t0_detected <- TRUE
  }
  if (t0_sec < 0) {
    warning("Time 0 (", round(t0_sec, 1),
            " s) is before the recording start; clamping to 0.")
    t0_sec <- 0
  }

  ## --- background, decay correction, calibration ----------------------------
  if (!is.null(bg_idx)) {
    # auto-detect: the baseline accumulated up to the rise onset
    background <- mean(coinc[bg_idx])
    n_background <- length(bg_idx)
  } else {
    # supplied pet_start: average everything before time 0
    bg_mask <- secs < t0_sec
    n_background <- sum(bg_mask)
    if (!any(bg_mask)) {
      warning("No samples before the supplied time 0; background set to 0.")
      background <- 0
    } else {
      background <- mean(coinc[bg_mask])
    }
  }

  tau  <- secs - t0_sec
  corr <- (coinc - background) * exp(lambda * tau) * calibration_factor

  ## --- cut before time 0 ----------------------------------------------------
  keep <- tau >= 0
  ord  <- order(tau[keep])
  tau_k  <- tau[keep][ord]
  corr_k <- corr[keep][ord]

  if (is.null(frame_scheme)) {
    ## --- no framing: return the corrected samples as-is ---------------------
    n <- length(tau_k)
    value <- corr_k
    if (zero_first_frame && n) value[1] <- 0
    # frame boundaries that tile the timeline (midpoints between samples)
    if (n >= 2) {
      mids        <- (tau_k[-1] + tau_k[-n]) / 2
      frame_start <- c(0, mids)
      frame_end   <- c(mids, tau_k[n] + (tau_k[n] - mids[n - 1L]))
    } else {
      frame_start <- rep(0, n)
      frame_end   <- tau_k
    }
    tac <- tibble::tibble(
      time        = tau_k,
      activity    = value,
      frame_start = frame_start,
      frame_end   = frame_end,
      frame_dur   = frame_end - frame_start
    )
  } else {
    ## --- build frame edges from the scheme ----------------------------------
    edges <- 0
    prev <- 0
    for (r in seq_len(nrow(frame_scheme))) {
      edges <- c(edges, seq(prev + frame_scheme$width[r], frame_scheme$end[r],
                            by = frame_scheme$width[r]))
      prev <- frame_scheme$end[r]
    }
    edges <- unique(edges)
    nfr <- length(edges) - 1L

    ## --- resample (average corrected samples per frame) ---------------------
    bin <- findInterval(tau_k, edges, rightmost.closed = FALSE)  # 1..nfr inside scheme
    inside <- bin >= 1L & bin <= nfr
    value <- tapply(corr_k[inside], factor(bin[inside], levels = seq_len(nfr)),
                    mean, simplify = TRUE)
    value <- as.numeric(value)  # NA where a frame had no samples

    # keep only frames that fall within the acquired data
    last_data_frame <- findInterval(max(tau_k), edges, rightmost.closed = TRUE)
    last_data_frame <- min(last_data_frame, nfr)
    idx <- seq_len(last_data_frame)

    if (zero_first_frame && length(idx)) value[idx[1]] <- 0

    tac <- tibble::tibble(
      time        = (edges[idx] + edges[idx + 1L]) / 2,
      activity    = value[idx],
      frame_start = edges[idx],
      frame_end   = edges[idx + 1L],
      frame_dur   = edges[idx + 1L] - edges[idx]
    )
  }

  ## --- assemble result ------------------------------------------------------
  list(
    tac                = tac,
    calibration_factor = calibration_factor,
    background         = background,
    isotope            = isotope,
    half_life          = half_life,
    lambda             = lambda,
    date               = as.Date(ts[1]),
    start_time         = format(ts[1], "%H:%M:%OS3"),
    acq_start          = ts[1],
    acq_end            = ts[length(ts)],
    acq_duration       = secs[length(secs)],
    pet_start          = ts[1] + t0_sec,
    t0_seconds         = t0_sec,
    t0_detected        = t0_detected,
    lead               = lead_used,
    n_raw              = length(secs),
    n_background       = n_background,
    file               = normalizePath(file),
    frame_scheme       = frame_scheme,
    raw = tibble::tibble(time = secs, coincidence = coinc,
                         singles1 = singles1, singles2 = singles2)
  )
}
