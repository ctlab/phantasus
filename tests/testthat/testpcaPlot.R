context("PCA Plot")

test_that("pcaPlot finishes with result", {
  load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
  expect_is(calcPCA(es), "json")
})

test_that("pcaPlot results in columnsXcolumns matrix", {
  load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
  expect_equal(nrow(jsonlite::fromJSON(calcPCA(es))$pca),
                length(sampleNames(es)))
})
