context("Create ExpressionSet")
library(Biobase)

test_that("createES finishes with result", {
  load(file = "testdata/GSE27112-GPL6103.rda")
  expect_is(createES(data = exprs(es), pData = pData(es),
                     fData = fData(es),
                     varLabels = varLabels(es),
                     fvarLabels = fvarLabels(es)), "ExpressionSet")
})
