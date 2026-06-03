# Each test corrects a simulated recording into a fresh temp output folder.

new_outdir <- function() {
  d <- tempfile("outputs_")
  dir.create(d)
  d
}

test_that("isotope is required", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  expect_error(swisstrace_process(path, calibration_factor = 0.4),
               "isotope.*required")
})

test_that("the corrected curve and QC plot are written", {
  path <- make_sim_crv()
  out  <- new_outdir()
  on.exit(unlink(c(path, out), recursive = TRUE))

  written <- suppressMessages(
    swisstrace_process(path, calibration_factor = 0.425, isotope = "F18",
                          output_folder = out))

  expect_true(file.exists(written$corrected_crv))
  expect_true(file.exists(written$plot_png))
  expect_match(written$corrected_crv, file.path("Corrected_PMOD", ".*_corrected\\.crv$"))
  expect_match(written$plot_png, file.path("Plots", ".*_corrected\\.png$"))

  # the written curve carries the expected header and matches res$tac
  hdr <- readLines(written$corrected_crv, n = 1)
  expect_match(hdr, "value\\[kBq/cc\\]")
  body <- read.table(written$corrected_crv, skip = 1, sep = "\t")
  res  <- swisstrace_correct(path, calibration_factor = 0.425, isotope = "F18")
  expect_equal(nrow(body), nrow(res$tac))
  expect_equal(body$V2, res$tac$activity, tolerance = 1e-6)
})

test_that("no BIDS output is written unless a subject is given", {
  path <- make_sim_crv()
  out  <- new_outdir()
  on.exit(unlink(c(path, out), recursive = TRUE))
  written <- suppressMessages(
    swisstrace_process(path, calibration_factor = 0.4, isotope = "F18",
                          output_folder = out))
  expect_null(written$bids_tsv)
  expect_false(dir.exists(file.path(out, "BIDS")))
})

test_that("ses without sub is an error", {
  path <- make_sim_crv()
  out  <- new_outdir()
  on.exit(unlink(c(path, out), recursive = TRUE))
  expect_error(
    swisstrace_process(path, calibration_factor = 0.4, isotope = "F18",
                          output_folder = out, ses = "01"),
    "without `sub`")
})

test_that("BIDS PET blood files are written under sub-/ses- with the right metadata", {
  path <- make_sim_crv()
  out  <- new_outdir()
  on.exit(unlink(c(path, out), recursive = TRUE))

  written <- suppressMessages(
    swisstrace_process(path, calibration_factor = 0.425, isotope = "F18",
                          output_folder = out, sub = "01", ses = "02"))

  expect_match(written$bids_tsv,
               file.path("BIDS", "sub-01", "ses-02", "pet",
                         "sub-01_ses-02_recording-autosampler_blood\\.tsv$"))
  expect_true(file.exists(written$bids_tsv))
  expect_true(file.exists(written$bids_json))

  tsv <- readLines(written$bids_tsv, n = 1)
  expect_equal(tsv, "time\twhole_blood_radioactivity")

  meta <- jsonlite::read_json(written$bids_json)
  expect_false(meta$PlasmaAvail)
  expect_true(meta$WholeBloodAvail)
  expect_false(meta$MetaboliteAvail)
  expect_false(meta$DispersionCorrected)
  expect_equal(meta$whole_blood_radioactivity$Units, "kBq/mL")
})

test_that("a subject with no session omits the ses- entity", {
  path <- make_sim_crv()
  out  <- new_outdir()
  on.exit(unlink(c(path, out), recursive = TRUE))
  written <- suppressMessages(
    swisstrace_process(path, calibration_factor = 0.4, isotope = "F18",
                          output_folder = out, sub = "07"))
  expect_match(written$bids_tsv,
               file.path("sub-07", "pet", "sub-07_recording-autosampler_blood\\.tsv$"))
  expect_false(grepl("ses-", written$bids_tsv))
})

test_that("output_folder defaults to the source file's directory", {
  out  <- new_outdir()
  path <- make_sim_crv(dir = out)
  on.exit(unlink(out, recursive = TRUE))
  written <- suppressMessages(
    swisstrace_process(path, calibration_factor = 0.4, isotope = "F18"))
  expect_true(startsWith(normalizePath(written$corrected_crv),
                         normalizePath(out)))
})

test_that("the returned path list is invisible", {
  path <- make_sim_crv()
  out  <- new_outdir()
  on.exit(unlink(c(path, out), recursive = TRUE))
  expect_invisible(suppressMessages(
    swisstrace_process(path, calibration_factor = 0.4, isotope = "F18",
                          output_folder = out)))
})

# --- bids_dir ----------------------------------------------------------------

test_that("bids_dir places blood files beside a matching _pet image, inheriting entities", {
  path <- make_sim_crv()
  bids <- new_outdir()
  on.exit(unlink(c(path, bids), recursive = TRUE))

  img_dir <- file.path(bids, "sub-01", "ses-02", "pet")
  dir.create(img_dir, recursive = TRUE)
  file.create(file.path(img_dir, "sub-01_ses-02_trc-PBR28_rec-acdyn_pet.nii.gz"))

  written <- suppressMessages(
    swisstrace_process(path, calibration_factor = 0.425, isotope = "F18",
                       output_folder = file.path(bids, "derivatives"),
                       sub = "01", ses = "02", bids_dir = bids))

  expect_equal(normalizePath(dirname(written$bids_tsv)), normalizePath(img_dir))
  expect_equal(basename(written$bids_tsv),
               "sub-01_ses-02_trc-PBR28_rec-acdyn_recording-autosampler_blood.tsv")
  expect_true(file.exists(written$bids_json))
})

test_that("bids_dir without a matching image falls back to sub-/ses-/pet under it", {
  path <- make_sim_crv()
  bids <- new_outdir()
  on.exit(unlink(c(path, bids), recursive = TRUE))

  written <- suppressMessages(
    swisstrace_process(path, calibration_factor = 0.4, isotope = "F18",
                       output_folder = file.path(bids, "derivatives"),
                       sub = "03", ses = "01", bids_dir = bids))

  expect_equal(normalizePath(dirname(written$bids_tsv)),
               normalizePath(file.path(bids, "sub-03", "ses-01", "pet")))
  expect_equal(basename(written$bids_tsv),
               "sub-03_ses-01_recording-autosampler_blood.tsv")
})

test_that("bids_dir requires sub", {
  path <- make_sim_crv()
  bids <- new_outdir()
  on.exit(unlink(c(path, bids), recursive = TRUE))
  expect_error(
    swisstrace_process(path, calibration_factor = 0.4, isotope = "F18",
                       bids_dir = bids),
    "without `sub`")
})

test_that("a sub label does not match a longer subject's image (sub-1 vs sub-10)", {
  path <- make_sim_crv()
  bids <- new_outdir()
  on.exit(unlink(c(path, bids), recursive = TRUE))

  d <- file.path(bids, "sub-10", "pet")
  dir.create(d, recursive = TRUE)
  file.create(file.path(d, "sub-10_pet.nii.gz"))

  written <- suppressMessages(
    swisstrace_process(path, calibration_factor = 0.4, isotope = "F18",
                       output_folder = file.path(bids, "derivatives"),
                       sub = "1", bids_dir = bids))

  # no sub-1 image exists, so it must fall back rather than land in sub-10/
  expect_equal(normalizePath(dirname(written$bids_tsv)),
               normalizePath(file.path(bids, "sub-1", "pet")))
})
