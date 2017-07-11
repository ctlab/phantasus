context("PCA Plot")

test_that("pcaPlot finishes with result", {
  load(file = "../../data/GSE27112-GPL6103.rda")
  expect_is(pcaPlot(es), "json")
})
