# Batch conversion driven by a simulated manifest. No real data is used.

batch_outdir <- function() {
  d <- tempfile("batch_")
  dir.create(d)
  d
}

write_manifest <- function(df, dir, ext = "csv") {
  path <- tempfile("manifest_", tmpdir = dir, fileext = paste0(".", ext))
  if (ext == "csv") {
    utils::write.csv(df, path, row.names = FALSE)
  } else {
    utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  path
}

cal_dates  <- c("2025-09-04", "2026-01-15", "2026-02-20")
cal_values <- c(0.198, 0.41, 0.425)

test_that("rows with a calibration_factor are converted without a lookup", {
  out <- batch_outdir()
  f1  <- make_sim_crv()
  f2  <- make_sim_crv(seed = 2)
  on.exit(unlink(c(out, f1, f2), recursive = TRUE))

  man <- data.frame(filename = c(f1, f2), isotope = "F18",
                    pet_start = c("11:16:30", ""),       # 2nd: blank -> auto-detect
                    calibration_factor = c(0.425, 0.40))
  res <- suppressMessages(swisstrace_convert_batch(write_manifest(man, out),
                                                   output_folder = out))

  expect_equal(res$status, c("ok", "ok"))
  expect_equal(res$calibration_factor, c(0.425, 0.40))
  expect_true(file.exists(file.path(out, "Corrected_PMOD",
              paste0(sub("\\.crv$", "", basename(f1)), "_corrected.crv"))))
  expect_false(file.exists(file.path(out, "calibration_factors.tsv")))   # no lookup
})

test_that("missing factors are looked up and a calibration log is written", {
  out <- batch_outdir()
  f1  <- make_sim_crv()
  on.exit(unlink(c(out, f1), recursive = TRUE))

  man <- data.frame(filename = f1, isotope = "F18", pet_start = "11:16:30")
  res <- suppressMessages(swisstrace_convert_batch(
    write_manifest(man, out), output_folder = out,
    cal_dates = cal_dates, cal_values = cal_values))

  expect_equal(res$status, "ok")
  expect_equal(res$calibration_factor, 0.425)         # study 2026-02-24 -> 2026-02-20

  log <- file.path(out, "calibration_factors.tsv")
  expect_true(file.exists(log))
  ll <- utils::read.delim(log)
  expect_named(ll, c("filename", "date", "calibration_factor", "gap_days"))
  expect_equal(ll$calibration_factor, 0.425)
  expect_equal(as.character(ll$date), "2026-02-24")
  expect_equal(ll$gap_days, 4)
})

test_that("a missing factor with no calibration log given is an error", {
  out <- batch_outdir()
  f1  <- make_sim_crv()
  on.exit(unlink(c(out, f1), recursive = TRUE))
  man <- data.frame(filename = f1, isotope = "F18")
  expect_error(
    suppressMessages(swisstrace_convert_batch(write_manifest(man, out), output_folder = out)),
    "cal_dates")
})

test_that("one bad row is reported as an error without aborting the others", {
  out <- batch_outdir()
  f1  <- make_sim_crv()
  on.exit(unlink(c(out, f1), recursive = TRUE))

  man <- data.frame(filename = c(f1, "no_such_file.crv"), isotope = "F18",
                    calibration_factor = c(0.425, 0.40))
  res <- suppressMessages(swisstrace_convert_batch(write_manifest(man, out),
                                                   output_folder = out))
  expect_equal(res$status, c("ok", "error"))
  expect_true(!is.na(res$message[2]))
})

test_that("sub/ses columns with bids_dir produce per-row BIDS output", {
  out  <- batch_outdir()
  bids <- batch_outdir()
  f1   <- make_sim_crv()
  on.exit(unlink(c(out, bids, f1), recursive = TRUE))

  man <- data.frame(filename = f1, isotope = "F18", pet_start = "11:16:30",
                    calibration_factor = 0.425, sub = "01", ses = "02")
  res <- suppressMessages(swisstrace_convert_batch(write_manifest(man, out),
                                                   output_folder = out, bids_dir = bids))
  expect_equal(res$status, "ok")
  expect_true(file.exists(file.path(bids, "sub-01", "ses-02", "pet",
                                    "sub-01_ses-02_recording-autosampler_blood.tsv")))
})

test_that("a .tsv manifest and a data.frame manifest both work", {
  out <- batch_outdir()
  f1  <- make_sim_crv()
  on.exit(unlink(c(out, f1), recursive = TRUE))
  man <- data.frame(filename = f1, isotope = "F18", calibration_factor = 0.425)

  from_tsv <- suppressMessages(swisstrace_convert_batch(
    write_manifest(man, out, ext = "tsv"), output_folder = out))
  expect_equal(from_tsv$status, "ok")

  from_df <- suppressMessages(swisstrace_convert_batch(man, output_folder = out))
  expect_equal(from_df$status, "ok")
})

test_that("a manifest without a filename column errors", {
  out <- batch_outdir()
  on.exit(unlink(out, recursive = TRUE))
  man <- data.frame(isotope = "F18", calibration_factor = 0.425)
  expect_error(
    suppressMessages(swisstrace_convert_batch(write_manifest(man, out), output_folder = out)),
    "missing required column")
})

test_that("an unsupported manifest extension errors", {
  out <- batch_outdir()
  p   <- tempfile(fileext = ".json")
  file.create(p)
  on.exit(unlink(c(out, p), recursive = TRUE))
  expect_error(swisstrace_convert_batch(p, output_folder = out), "Unsupported manifest")
})
