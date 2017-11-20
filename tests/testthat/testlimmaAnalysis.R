context("Limma Analysis")

test_that("limmaAnalysis finishes with result", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(limmaAnalysis(es,
                          fieldValues = c(rep("A", 3), rep("B", 2))), "json")
})

test_that("limmaAnalysis works when there is only one phenotype attribute", {
    load(file = "testdata/GSE27112-GPL6103.rda")
    pData(es) <- pData(es)[, "title", drop=F]
    expect_is(limmaAnalysis(es,
                            fieldValues = c(rep("A", 3), rep("B", 2))), "json")
})
