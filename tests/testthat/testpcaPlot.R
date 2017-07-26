context("PCA Plot")

test_that("pcaPlot finishes with result", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(pcaPlot(es), "json")
})

test_that("pcaPlot results in columnsXcolumns matrix", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(nrow(jsonlite::fromJSON(pcaPlot(es))$pca), length(sampleNames(es)))
})
