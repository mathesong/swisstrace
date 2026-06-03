test_that("a simulated raw recording passes the assertion", {
  path <- make_sim_crv()
  on.exit(unlink(path))
  expect_true(assert_raw_crv(path))
  expect_invisible(assert_raw_crv(path))
})

test_that("a corrected/processed file is rejected", {
  path <- tempfile(fileext = ".crv")
  on.exit(unlink(path))
  write_corrected_crv(path)
  expect_error(assert_raw_crv(path), "corrected/processed")
})

test_that("a non-numeric first line is rejected", {
  path <- tempfile(fileext = ".crv")
  on.exit(unlink(path))
  writeLines(c("patient header line", "2026 2 24 11 15 21.256\t33\t506\t944"), path)
  expect_error(assert_raw_crv(path), "looks like a corrected/processed")
})

test_that("an empty file is rejected", {
  path <- tempfile(fileext = ".crv")
  on.exit(unlink(path))
  writeLines(character(), path)
  expect_error(assert_raw_crv(path), "empty")
})

test_that("a missing file errors", {
  expect_error(assert_raw_crv(tempfile(fileext = ".crv")))
})
