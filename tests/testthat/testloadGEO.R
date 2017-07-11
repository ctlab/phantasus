context("Load GEO")

test_that("loadGEO works", {
  expect_is(loadGEO("GSE27112"), "json")
  expect_is(loadGEO("GSE27112-GPL6885"), "json")
  expect_is(loadGEO("GSE14308"), "json")
  expect_is(loadGEO("GDS4885"), "json")
})

