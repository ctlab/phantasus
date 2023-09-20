context("GCT")

test_that("write.gct and read.gct work", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
    gctFile <- tempfile(fileext = ".gct")
    dir.create(dirname(gctFile), recursive = TRUE)
    write.gct(es, gctFile)
    es2 <- read.gct(gctFile)
    expect_equal(dim(es2), dim(es))
})

test_that("write.gct and read.gct work for gzip files", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
    gctFile <- tempfile(fileext = ".gct.gz")
    dir.create(dirname(gctFile), recursive = TRUE)
    write.gct(es, gctFile, gzip = T)
    es2 <- read.gct(gctFile)
    expect_equal(dim(es2), dim(es))
})

test_that("read.gct works (simple version)", {
    es <- read.gct(system.file("testdata/test.gct", package="phantasus"))
    expect_equal(dim(exprs(es)), c(5, 7))
    expect_equal(colnames(es)[1], "s1")
    expect_equal(rownames(es)[1], "1415670_at")
    expect_true(all(c("Gene symbol", "Gene ID") %in% colnames(fData(es))), setdiff(c("Gene symbol", "Gene ID"), colnames(fData(es))))
    expect_true(all(c("condition") %in% colnames(pData(es))), setdiff(c("condition"), colnames(pData(es))))
})
