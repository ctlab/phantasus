context("Util functions")

test_that("prepareData replaces all NAs and filters fully NA rows", {
  data <- prepareData(read.gct("testdata/centers.gct"))
  expect_equal(nrow(data), 5)
})

test_that("data preparation successfully eliminates na-rows when at first rows are specified", {
  load("testdata/aa.rda")
  data <- prepareData(aa, rows = unlist(jsonlite::read_json("testdata/rows")))
  expect_equal(nrow(data), 10)
})
