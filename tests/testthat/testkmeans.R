context("Kmeans")

test_that("kmeans finishes with result", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(phantasus::kmeans(es, k = 10), "json")
})
