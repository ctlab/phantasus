context("Util functions")

test_that("prepareData replaces all NAs and filters fully NA rows", {
  data <- prepareData(read.gct(system.file("testdata/centers.gct", package="phantasus")))
  expect_equal(nrow(data), 5)
})

test_that("data preparation successfully eliminates na-rows when at first rows are specified", {
  load(system.file("testdata/aa.rda", package="phantasus"))
  data <- prepareData(aa, rows = unlist(jsonlite::read_json(system.file("testdata/rows", package="phantasus"))))
  expect_equal(nrow(data), 10)
})
