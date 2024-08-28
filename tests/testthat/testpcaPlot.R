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

test_that("pcaPlot drops NA rows", {
    es <- readGct(system.file("testdata/centers.gct", package="phantasus"))
    pca <- jsonlite::fromJSON(calcPCA(es))
    expect_equal(nrow(exprs(es)), 6)
    expect_equal(nrow(pca$pca), 5)
})
