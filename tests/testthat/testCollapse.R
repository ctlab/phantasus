context('Collapse dataset')

test_that("Collapse dataset selectOne works correctly", {
    es <- read.gct(system.file("testdata/collapse_dataset_one.gct", package="phantasus"))
    newEs <- collapseDatasetImpl(es, selectOne = TRUE, fn = median, fields = c('Gene ID'), removeEmpty = FALSE)
    expect_equal(dim(exprs(newEs)), c(3,7))
    expect_equal(unname(dim(newEs)), dim(exprs(newEs)))
    expect_equal(fData(newEs)[1,][['Gene symbol']], 'Copg1')
    expect_true(all(isUnique(fData(newEs)[['Gene ID']])))
    expect_true(all(exprs(es)[5,] == exprs(newEs)[3,]))
    expect_true(all(exprs(es)[4,] == exprs(newEs)[2,]))
    expect_true(all(exprs(es)[1,] == exprs(newEs)[1,]))
})


test_that("Collapse dataset selectOne all unique", {
    es <- read.gct(system.file("testdata/collapse_dataset_one.gct", package="phantasus"))
    newEs <- collapseDatasetImpl(es, selectOne = TRUE, fn = median, fields = c('id','Gene ID', 'Gene symbol'), removeEmpty = FALSE)
    expect_equal(dim(exprs(newEs)), c(5,7))
    expect_equal(unname(dim(newEs)), dim(exprs(newEs)))
    expect_true(all(exprs(es) == exprs(newEs)))
})

test_that("Collapse dataset by rows", {
    es <- read.gct(system.file("testdata/collapse_dataset_one.gct", package="phantasus"))
    newEs <- collapseDatasetImpl(es, fn = median, fields = c('Gene ID'), removeEmpty = FALSE)
    expect_equal(ncol(fData(newEs)), 1)
    expect_equal(nrow(fData(newEs)), 3)
    expect_equal(dim(exprs(newEs)), c(3,7))
    expect_true(all(isUnique(fData(newEs)[['Gene ID']])))
    expect_true(all(exprs(es)[5,] == exprs(newEs)[3,]))
    expect_true(all(exprs(es)[4,] == exprs(newEs)[2,]))
    expect_equal(apply(exprs(es)[c(1,2,3),], 2, median), exprs(newEs)[1,])


    newEsMultiple <- collapseDatasetImpl(es, fn = mean, fields = c('Gene ID', 'Gene symbol'), removeEmpty = FALSE)
    expect_equal(ncol(fData(newEsMultiple)), 2)
    expect_equal(apply(exprs(es)[c(1,2,3),], 2, mean), exprs(newEsMultiple)[1,])
    expect_true(all(isUnique(paste(fData(newEsMultiple)[['Gene ID']], fData(newEsMultiple)[['Gene symbol']], sep="//r"))))
})

test_that("Collapse dataset colnames doesn't reoder", {
    es <- read.gct(system.file("testdata/collapse_dataset_one.gct", package="phantasus"))
    newEs <- collapseDatasetImpl(es, fn = median, fields = c('Gene ID', 'id', 'Gene symbol'), removeEmpty = FALSE)
    expect_equal(colnames(fData(newEs)), colnames(fData(es)))

    newEs <- collapseDatasetImpl(es, fn = median, fields = c('Gene ID', 'Gene symbol'), removeEmpty = FALSE)
    expect_lt(match('Gene symbol', colnames(fData(newEs))), match('Gene ID', colnames(fData(newEs))))
})


test_that("Collapse dataset by columns", {
    es <- read.gct(system.file("testdata/collapse_dataset_one.gct", package="phantasus"))
    newEs <- collapseDatasetImpl(es, isRows = FALSE, fn = median, fields = c('condition'), removeEmpty = FALSE)
    expect_equal(ncol(pData(newEs)), 1)
    expect_equal(nrow(pData(newEs)), 2)
    expect_equal(dim(exprs(newEs)), c(5,2))
    expect_true(all(isUnique(pData(newEs)[['condition']])))
    expect_equal(apply(exprs(es)[,c(1,2,3,4)], 1, median), exprs(newEs)[,1])
    expect_equal(apply(exprs(es)[,c(5,6,7)], 1, median), exprs(newEs)[,2])
})

test_that("Collapse strip unannotated rows", {
    old_conf <- Sys.getenv("R_CONFIG_ACTIVE")
    Sys.setenv(R_CONFIG_ACTIVE = "test_real_geo")
    ess <- getGSE('GSE53986')[[1]]
    unstripped <- collapseDatasetImpl(ess, selectOne = TRUE, fn = median, fields = c('Gene ID'), removeEmpty = FALSE)
    stripped <- collapseDatasetImpl(ess, selectOne = TRUE, fn = median, fields = c('Gene ID'), removeEmpty = TRUE)
    expect_equal(nrow(stripped) + 1, nrow(unstripped))
    Sys.setenv(R_CONFIG_ACTIVE = old_conf)
})
