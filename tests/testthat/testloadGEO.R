context("Load GEO and its utils")
library(jsonlite)
library(Biobase)

test_that("loadGEO finishes with result", {
  expect_is(loadGEO("GSE27112"), "json")
  expect_is(loadGEO("GSE27112-GPL6885"), "json")
  expect_is(loadGEO("GSE14308"), "json")
  expect_is(loadGEO("GDS4885"), "json")
})

test_that("checkGPLs counts gpls correctly", {
  expect_equal(fromJSON(checkGPLs("GSE14308")), c("GSE14308"))
  expect_equal(fromJSON(checkGPLs("GDS4885")), c("GDS4885"))
  expect_length(fromJSON(checkGPLs("GSE27112")), 2)
  expect_length(fromJSON(checkGPLs("GSE10000")), 2)
  expect_length(fromJSON(checkGPLs("GDS101")), 0)
  expect_length(fromJSON(checkGPLs("GSE201")), 0)
})
