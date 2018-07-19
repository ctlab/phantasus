context("GCT")

test_that("write.gct and read.gct work", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
    gctFile <- tempfile(fileext = ".gct")
    write.gct(es, gctFile)
    es2 <- read.gct(gctFile)
    expect_equal(dim(es2), dim(es))
})

test_that("write.gct and read.gct work for gzip files", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
    gctFile <- tempfile(fileext = ".gct.gz")
    write.gct(es, gctFile, gzip = T)
    es2 <- read.gct(gctFile)
    expect_equal(dim(es2), dim(es))
})
