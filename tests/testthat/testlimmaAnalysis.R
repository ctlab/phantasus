context("Limma Analysis")

test_that("limmaAnalysis finishes with result", {
  load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
  expect_is(limmaAnalysis(es,
                          fieldValues = c(rep("A", 3), rep("B", 2))), "json")
})

test_that("limmaAnalysis works when there is only one phenotype attribute", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
    pData(es) <- pData(es)[, "title", drop=F]
    expect_is(limmaAnalysis(es,
                            fieldValues = c(rep("A", 3), rep("B", 2))), "json")
})


test_that("limmaAnalysisImpl works", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
    de <- limmaAnalysisImpl(es, fieldValues = c(rep("A", 2), rep("B", 2), NA))

    expect_equal(de$logFC[1], mean(exprs(es)[1, 3:4]) - mean(exprs(es)[1, 1:2]))
})

# Limma works for full dataset only
#test_that("limmaAnalysisImpl works for subsamples", {
#    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
#    de1 <- limmaAnalysisImpl(es, rows=seq_len(nrow(es)), columns=seq_len(ncol(es)),
#                            fieldValues = c(rep("A", 2), NA, rep("B", 2)))
#
#    de2 <- limmaAnalysisImpl(es, rows=seq_len(nrow(es)), columns=c(1:2, 4:5),
#                             fieldValues = c(rep("A", 2), rep("B", 2)))
#    expect_equal(de1$t, de2$t)
#})
