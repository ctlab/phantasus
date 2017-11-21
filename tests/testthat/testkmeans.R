context("Kmeans")
library(jsonlite)

test_that("kmeans finishes with result", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(performKmeans(es, k = 10), "json")
})

test_that("kmeans works with datasets with lots of NAs", {
  es <- read.gct("testdata/centers.gct")
  expect_length(fromJSON(performKmeans(es, k = 2)), nrow(es))
})
