context("Limma Analysis")

test_that("limmaAnalysis finishes with result", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(limmaAnalysis(es,
                          fieldValues = c(rep("A", 3), rep("B", 2))), "json")
})
