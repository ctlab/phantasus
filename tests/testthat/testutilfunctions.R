context("Util functions")
library(rUtils)

test_that("prepareData replaces all NAs and filters fully NA rows", {
  data <- prepareData(es)
  expect_equal(nrow(data), 5)
})
