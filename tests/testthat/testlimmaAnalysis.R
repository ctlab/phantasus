context("Limma Analysis")

test_that("limmaAnalysis finishes with result", {
  load(file = "../../data/GSE27112-GPL6103.rda")
  expect_is(limmaAnalysis(es, fieldValues = c(rep("A", 5), rep("B", 5))), "json")
})
