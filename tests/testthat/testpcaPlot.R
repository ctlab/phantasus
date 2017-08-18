context("PCA Plot")

test_that("pcaPlot finishes with result", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(calcPCA(es), "json")
})

test_that("pcaPlot results in columnsXcolumns matrix", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_equal(nrow(jsonlite::fromJSON(calcPCA(es))$pca), length(sampleNames(es)))
})
