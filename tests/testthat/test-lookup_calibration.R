cal_dates  <- as.Date(c("2025-09-04", "2026-01-15", "2026-02-20"))
cal_values <- c(0.198, 0.41, 0.425)

test_that("swisstrace_date reads the acquisition date from a raw file", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  expect_equal(swisstrace_date(path), as.Date("2026-02-24"))
})

test_that("method 'last' picks the most recent calibration on or before the study", {
  out <- lookup_calibration(date = "2026-02-24", cal_dates, cal_values, method = "last")
  expect_equal(out$calibration_factor, 0.425)
  expect_equal(out$cal_date, as.Date("2026-02-20"))
  expect_equal(out$gap_days, 4)            # study minus calibration date
})

test_that("method 'last' uses the calibration on the same day if present", {
  out <- lookup_calibration(date = "2026-01-15", cal_dates, cal_values, method = "last")
  expect_equal(out$calibration_factor, 0.41)
  expect_equal(out$gap_days, 0)
})

test_that("method 'nearest' can match a later calibration", {
  out <- lookup_calibration(date = "2026-01-10", cal_dates, cal_values, method = "nearest")
  expect_equal(out$calibration_factor, 0.41)   # 2026-01-15 is closest
  expect_equal(out$gap_days, -5)
})

test_that("method 'exact' only matches an exact date", {
  hit  <- lookup_calibration(date = "2026-01-15", cal_dates, cal_values, method = "exact")
  expect_equal(hit$calibration_factor, 0.41)
  expect_warning(
    miss <- lookup_calibration(date = "2026-01-16", cal_dates, cal_values, method = "exact"),
    "no calibration factor")
  expect_true(is.na(miss$calibration_factor))
})

test_that("a study before any calibration yields NA with a warning under 'last'", {
  expect_warning(
    out <- lookup_calibration(date = "2025-01-01", cal_dates, cal_values, method = "last"),
    "on or before")
  expect_true(is.na(out$calibration_factor))
})

test_that("max_gap blanks out matches that are too far away", {
  expect_warning(
    out <- lookup_calibration(date = "2026-06-01", cal_dates, cal_values,
                              method = "last", max_gap = 30),
    "within range")
  expect_true(is.na(out$calibration_factor))
})

test_that("lookup is vectorised over study dates", {
  out <- lookup_calibration(date = c("2025-09-10", "2026-02-24"), cal_dates, cal_values)
  expect_equal(nrow(out), 2)
  expect_equal(out$calibration_factor, c(0.198, 0.425))
})

test_that("a .crv path (filename) and a swisstrace_correct() result both resolve to a date", {
  path <- make_sim_crv()                     # acquisition date 2026-02-24
  on.exit(unlink(path))
  from_path <- lookup_calibration(filename = path, cal_dates, cal_values)
  expect_equal(from_path$calibration_factor, 0.425)

  res <- swisstrace_correct(path, calibration_factor = 0.425, isotope = "F18")
  from_res <- lookup_calibration(date = res, cal_dates, cal_values)
  expect_equal(from_res$study_date, as.Date("2026-02-24"))
  expect_equal(from_res$calibration_factor, 0.425)
})

test_that("a YYYY/MM/DD date string is accepted like YYYY-MM-DD", {
  slash <- lookup_calibration(date = "2026/02/24", cal_dates, cal_values)
  dash  <- lookup_calibration(date = "2026-02-24", cal_dates, cal_values)
  expect_equal(slash$study_date, as.Date("2026-02-24"))
  expect_equal(slash$calibration_factor, dash$calibration_factor)
})

test_that("supplying neither filename nor date errors", {
  expect_error(
    lookup_calibration(cal_dates = cal_dates, cal_values = cal_values),
    "Supply a study")
})

test_that("supplying both filename and date errors", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  expect_error(
    lookup_calibration(filename = path, date = "2026-02-24", cal_dates, cal_values),
    "only one")
})

test_that("a non-existent .crv filename errors clearly", {
  expect_error(
    lookup_calibration(filename = "ABC_P009_D1.crv", cal_dates, cal_values),
    "File not found")
})

test_that("an unparseable date string errors clearly", {
  expect_error(
    lookup_calibration(date = "not_a_date", cal_dates, cal_values),
    "Could not interpret")
})

test_that("mismatched cal_dates / cal_values lengths error", {
  expect_error(
    lookup_calibration(date = "2026-02-24", cal_dates, cal_values[1:2]),
    "same length")
})
