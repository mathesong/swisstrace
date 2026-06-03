test_that("swisstrace_qc draws without error and returns its input invisibly", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  res <- swisstrace_correct(path, calibration_factor = 0.425, isotope = "F18")

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp)
  on.exit(grDevices::dev.off(), add = TRUE)

  out <- swisstrace_qc(res)
  expect_identical(out, res)
  expect_invisible(swisstrace_qc(res))
  expect_true(file.exists(tmp))
})
