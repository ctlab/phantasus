context("Kmeans")
library(jsonlite)

test_that("kmeans finishes with result", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(kmeans(es, k = 10), "json")
})

test_that("kmeans works with datasets with lots of NAs, so there is need to filter rows", {
  es <- read.gct("testdata/centers.gct")
  expect_length(fromJSON(kmeans(es, k = 2)), 5)
})
