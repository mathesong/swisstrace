# --- input validation --------------------------------------------------------

test_that("isotope is required when no half_life is given", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  expect_error(swisstrace_correct(path, calibration_factor = 0.4),
               "isotope.*required")
})

test_that("an unknown isotope errors and lists the known ones", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  expect_error(swisstrace_correct(path, calibration_factor = 0.4, isotope = "Xx99"),
               "Unknown isotope")
})

test_that("half_life can stand in for isotope", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  res <- swisstrace_correct(path, calibration_factor = 0.4, half_life = 6586.2)
  expect_equal(res$half_life, 6586.2)
  expect_true(is.na(res$isotope))
  expect_equal(res$lambda, log(2) / 6586.2)
})

test_that("isotope names are matched case/punctuation-insensitively", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  a <- swisstrace_correct(path, calibration_factor = 0.4, isotope = "F18")
  b <- swisstrace_correct(path, calibration_factor = 0.4, isotope = "f-18")
  expect_equal(a$half_life, b$half_life)
  expect_equal(a$half_life, 6586.2)
})

test_that("calibration_factor must be a single number", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  expect_error(swisstrace_correct(path, isotope = "F18"),
               "calibration_factor")
  expect_error(swisstrace_correct(path, calibration_factor = c(1, 2), isotope = "F18"),
               "single number")
})

test_that("a corrected .crv is rejected up front", {
  path <- tempfile(fileext = ".crv")
  on.exit(unlink(path))
  write_corrected_crv(path)
  expect_error(swisstrace_correct(path, calibration_factor = 0.4, isotope = "F18"),
               "corrected/processed")
})

# --- the correction arithmetic ----------------------------------------------

test_that("a constant decay-corrected curve yields constant calibrated activity", {
  cal <- 0.5
  g   <- make_exact_crv(background = 30, amp = 100, t0 = 40, half_life = 6586.2)
  on.exit(unlink(g$path))

  res <- swisstrace_correct(g$path, calibration_factor = cal,
                             half_life = g$half_life,
                             pet_start = g$t0, zero_first_frame = TRUE)

  # background recovered exactly, decay constant correct
  expect_equal(res$background, g$background)
  expect_equal(res$lambda, log(2) / g$half_life)
  expect_equal(res$t0_seconds, g$t0)
  expect_false(res$t0_detected)

  # every frame after the (zeroed) first equals amp * calibration
  act <- res$tac$activity
  expect_equal(act[1], 0)
  expect_equal(act[-1], rep(g$amp * cal, length(act) - 1L), tolerance = 1e-6)
})

test_that("pet_start works identically as POSIXct clock time or seconds", {
  g <- make_exact_crv(t0 = 40)
  on.exit(unlink(g$path))
  by_secs  <- swisstrace_correct(g$path, calibration_factor = 0.5,
                                  half_life = g$half_life, pet_start = g$t0)
  by_clock <- swisstrace_correct(g$path, calibration_factor = 0.5,
                                  half_life = g$half_life,
                                  pet_start = g$start + g$t0)
  expect_equal(by_secs$t0_seconds, by_clock$t0_seconds)
  expect_equal(by_secs$tac$activity, by_clock$tac$activity)
})

test_that("pet_start accepts an HH:MM:SS clock time on the recording's date", {
  g <- make_exact_crv(t0 = 40)   # recording starts 2026-01-01 00:00:00, so t0 = 00:00:40
  on.exit(unlink(g$path))
  by_str  <- swisstrace_correct(g$path, calibration_factor = 0.5,
                                half_life = g$half_life, pet_start = "00:00:40")
  by_secs <- swisstrace_correct(g$path, calibration_factor = 0.5,
                                half_life = g$half_life, pet_start = g$t0)
  expect_false(by_str$t0_detected)
  expect_equal(by_str$t0_seconds, 40)
  expect_equal(by_str$tac$activity, by_secs$tac$activity)
})

test_that("an unparseable pet_start clock time errors helpfully", {
  g <- make_exact_crv(t0 = 40)
  on.exit(unlink(g$path))
  expect_error(
    swisstrace_correct(g$path, calibration_factor = 0.5, half_life = g$half_life,
                       pet_start = "not a time"),
    "valid clock time")
})

test_that("background is the mean of pre-t0 samples for an explicit pet_start", {
  g <- make_exact_crv(background = 30, t0 = 40)
  on.exit(unlink(g$path))
  res <- swisstrace_correct(g$path, calibration_factor = 1,
                             half_life = g$half_life, pet_start = g$t0)
  expect_equal(res$n_background, 40L)   # samples at tau 0..39
  expect_equal(res$background, 30)
})

test_that("zero_first_frame toggles the first frame", {
  g <- make_exact_crv()
  on.exit(unlink(g$path))
  on  <- swisstrace_correct(g$path, calibration_factor = 0.5,
                             half_life = g$half_life, pet_start = g$t0,
                             zero_first_frame = TRUE)
  off <- swisstrace_correct(g$path, calibration_factor = 0.5,
                             half_life = g$half_life, pet_start = g$t0,
                             zero_first_frame = FALSE)
  expect_equal(on$tac$activity[1], 0)
  expect_gt(off$tac$activity[1], 0)
})

# --- output structure --------------------------------------------------------

test_that("the result has the documented shape", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  res <- swisstrace_correct(path, calibration_factor = 0.425, isotope = "F18")

  expect_s3_class(res$tac, "tbl_df")
  expect_named(res$tac, c("time", "activity", "frame_start", "frame_end", "frame_dur"))
  expect_equal(res$tac$frame_dur, res$tac$frame_end - res$tac$frame_start)
  expect_equal(res$tac$time, (res$tac$frame_start + res$tac$frame_end) / 2)

  expect_s3_class(res$raw, "tbl_df")
  expect_named(res$raw, c("time", "coincidence", "singles1", "singles2"))
  expect_equal(res$isotope, "F18")
  expect_s3_class(res$date, "Date")
  expect_equal(res$n_raw, nrow(res$raw))
})

test_that("a custom frame scheme drives the frame edges", {
  g <- make_exact_crv(total = 300)
  on.exit(unlink(g$path))
  res <- swisstrace_correct(g$path, calibration_factor = 0.5, half_life = g$half_life,
                             pet_start = g$t0,
                             frame_scheme = data.frame(width = 5, end = 100))
  expect_true(all(res$tac$frame_dur == 5))
  expect_lte(max(res$tac$frame_end), 100)
})

# --- automatic time-0 detection ---------------------------------------------

test_that("the bolus onset is detected and time 0 is placed a lead before it", {
  onset <- 60
  path  <- make_sim_crv(onset = onset, noise_sd = 2, seed = 7)
  on.exit(unlink(path))
  res <- swisstrace_correct(path, calibration_factor = 0.425, isotope = "F18")

  expect_true(res$t0_detected)
  expect_equal(res$lead, 20)
  # detected onset ~ 60 s, so time 0 ~ 40 s; allow a few samples of slack
  expect_equal(res$t0_seconds, onset - 20, tolerance = 4)
  expect_lt(abs(res$background - 32), 4)   # recovered baseline near the true 32
})

test_that("an insufficient pre-rise baseline shrinks the lead with a message", {
  onset <- 14
  path  <- make_sim_crv(onset = onset, total = 600, noise_sd = 2, seed = 3)
  on.exit(unlink(path))
  expect_message(
    res <- suppressWarnings(
      swisstrace_correct(path, calibration_factor = 0.4, isotope = "F18")),
    "lead"
  )
  expect_lt(res$lead, 20)
})
