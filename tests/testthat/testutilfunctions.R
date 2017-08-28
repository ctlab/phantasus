context("Util functions")

test_that("prepareData replaces all NAs and filters fully NA rows", {
  data <- prepareData(read.gct("testdata/centers.gct"))
  expect_equal(nrow(data), 5)
})
