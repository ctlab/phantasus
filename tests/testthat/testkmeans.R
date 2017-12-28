context("Kmeans")
library(jsonlite)

test_that("kmeans finishes with result", {
  load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
  expect_is(performKmeans(es, k = 10), "json")
})

test_that("kmeans works with datasets with lots of NAs", {
  es <- read.gct(system.file("testdata/centers.gct", package="phantasus"))
  expect_length(fromJSON(performKmeans(es, k = 2)), nrow(es))
})

test_that("kmeans returns only vector for a subset", {
    es <- read.gct(system.file("testdata/centers.gct", package="phantasus"))
    rows <- c(0, 2, 3, 5)
    expect_length(fromJSON(performKmeans(es, rows=rows, k = 2)), length(rows))
})
