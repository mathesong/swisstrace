# Helpers to fabricate raw twilite .crv data for the tests.
#
# No real recordings are bundled or read; everything here is simulated to mimic the
# shape of a genuine twilite file (tab-separated rows of
# "YYYY M D H M S.sss  coincidence  singles1  singles2", ~1 Hz, a flat baseline
# followed by a bolus rise and decay). The first column (coincidences) is the only
# channel used by the correction; the two singles channels are filler.

# Write raw rows to a .crv file from numeric vectors. `coinc`, `singles1`, `singles2`
# are per-sample counts; `start` is the POSIXct timestamp of the first sample, `step`
# the spacing in seconds, `frac` the fixed fractional-second offset real files carry.
write_raw_crv <- function(path, coinc,
                          start    = as.POSIXct("2026-02-24 11:15:21", tz = "UTC"),
                          singles1 = NULL,
                          singles2 = NULL,
                          step = 1, frac = 0.256,
                          round_counts = TRUE) {
  n <- length(coinc)
  if (is.null(singles1)) singles1 <- rep(550L, n)
  if (is.null(singles2)) singles2 <- rep(1000L, n)
  t  <- start + (seq_len(n) - 1L) * step
  lt <- as.POSIXlt(t, tz = "UTC")
  stamp <- sprintf("%d %d %d %d %d %.3f",
                   lt$year + 1900L, lt$mon + 1L, lt$mday, lt$hour, lt$min,
                   lt$sec + frac)
  # Genuine recordings store integer counts; the exact-arithmetic tests instead
  # write full-precision values so the round trip is lossless.
  cfmt <- if (round_counts) {
    function(x) format(as.integer(round(x)))
  } else {
    function(x) formatC(x, digits = 12, format = "g")
  }
  rows <- paste(stamp, cfmt(coinc),
                format(as.integer(round(singles1))),
                format(as.integer(round(singles2))), sep = "\t")
  writeLines(rows, path)
  invisible(path)
}

# Build a realistic-looking coincidence trace: flat baseline, a sharp bolus rise at
# `onset` s, then a washout combined with physical decay. Returns the per-second
# coincidence vector; pair it with write_raw_crv(). Deterministic given `seed`.
sim_coincidence <- function(secs, onset,
                            background = 32,
                            peak       = 220,
                            half_life  = 6586.2,
                            rise_tau   = 2,
                            washout    = 1200,
                            noise_sd   = 3,
                            seed       = 1) {
  set.seed(seed)
  lambda <- log(2) / half_life
  out  <- rep(background, length(secs))
  post <- secs >= onset
  dt   <- secs[post] - onset
  shape <- (1 - exp(-dt / rise_tau)) * exp(-dt / washout)  # rise then washout
  if (max(shape) > 0) shape <- shape / max(shape)
  out[post] <- background + (peak - background) * shape * exp(-lambda * dt)
  pmax(out + stats::rnorm(length(secs), 0, noise_sd), 0)
}

# Convenience: write a standard, well-behaved simulated recording to a temp .crv and
# return its path (long flat baseline so auto-detection has plenty to work with).
make_sim_crv <- function(onset = 60, total = 900, ..., dir = tempdir()) {
  secs  <- 0:total
  coinc <- sim_coincidence(secs, onset = onset, ...)
  path  <- tempfile(pattern = "twilite_", tmpdir = dir, fileext = ".crv")
  write_raw_crv(path, coinc)
  path
}

# A synthetic curve whose decay-CORRECTED signal is exactly constant after t0, so the
# end-to-end correction must return a flat activity of `amp * calibration`. Lets the
# tests check the arithmetic to numerical precision (no noise). Returns the .crv path
# plus the ground-truth pieces used to build it.
make_exact_crv <- function(background = 30, amp = 100, t0 = 40, total = 300,
                           half_life = 6586.2,
                           start = as.POSIXct("2026-01-01 00:00:00", tz = "UTC"),
                           dir = tempdir()) {
  lambda <- log(2) / half_life
  secs   <- 0:total
  coinc  <- rep(background, length(secs))
  post   <- secs >= t0
  coinc[post] <- background + amp * exp(-lambda * (secs[post] - t0))
  path <- tempfile(pattern = "exact_", tmpdir = dir, fileext = ".crv")
  # frac = 0 keeps the recording origin on `start` exactly (so a clock-time
  # pet_start lines up), and full-precision counts make the arithmetic lossless.
  write_raw_crv(path, coinc, start = start, frac = 0, round_counts = FALSE)
  list(path = path, background = background, amp = amp, t0 = t0,
       half_life = half_life, lambda = lambda, start = start, n = length(secs))
}

# Write a *corrected* (non-raw) .crv, to check that assertions reject it.
write_corrected_crv <- function(path) {
  writeLines(c("Corrected_&_calibrated_[kBq/cc]_>___time[seconds]\tvalue[kBq/cc]",
               "0.5\t0.0",
               "1.5\t2.30"),
             path)
  invisible(path)
}
