context("Load GEO and its utils")

test_that("loadGEO finishes with result", {
  expect_is(loadGEO("GSE27112"), "json")
  expect_is(loadGEO("GSE27112-GPL6885"), "json")
  expect_is(loadGEO("GSE14308"), "json")
  expect_is(loadGEO("GDS4885"), "json")
})

test_that("checkGPLs counts gpls correctly", {
  expect_equal(jsonlite::fromJSON(checkGPLs("GSE14308")), c("GSE14308"))
  expect_equal(length(jsonlite::fromJSON(checkGPLs("GSE27112"))), 2)
  expect_equal(length(jsonlite::fromJSON(checkGPLs("GSE10000"))), 2)
  expect_equal(jsonlite::fromJSON(checkGPLs("GDS4885")), c("GDS4885"))
  expect_equal(length(jsonlite::fromJSON(checkGPLs("GDS101"))), 0)
  expect_equal(length(jsonlite::fromJSON(checkGPLs("GSE201"))), 0)
})

