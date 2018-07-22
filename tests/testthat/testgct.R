context("GCT")

test_that("write.gct and read.gct work", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
    gctFile <- tempfile(fileext = ".gct")
    write.gct(es, gctFile)
    es2 <- read.gct(gctFile)
    expect_equal(dim(es2), dim(es))

    expect_true(all(colnames(fData(es)) %in% colnames(fData(es2))),
                info=setdiff( colnames(fData(es)), colnames(fData(es2))))

    expect_true(all(colnames(pData(es)) %in% colnames(pData(es2))),
                info=setdiff( colnames(pData(es)), colnames(pData(es2))))
})

test_that("write.gct and read.gct work for gzip files", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))
    gctFile <- tempfile(fileext = ".gct.gz")
    write.gct(es, gctFile, gzip = T)
    es2 <- read.gct(gctFile)
    expect_equal(dim(es2), dim(es))
})

test_that("read.gct works with duplicate row names", {
    expect_warning(es <- read.gct(system.file("testdata/dupl.gct", package="phantasus")))
})
