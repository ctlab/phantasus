context("Load GEO and its utils")
library(jsonlite)
library(Biobase)

test_that("loadGEO finishes with result", {
    options(phantasusMirrorPath = "https://genome.ifmo.ru/files/software/phantasus")

    x <- loadGEO("GSE27112")
    expect_is(x, "json")

    ess <- protolite::unserialize_pb(readBin(fromJSON(x), what="raw", n=100000000))

    expect_equal(length(ess), 2)

    x <- loadGEO("GSE27112-GPL6885")
    expect_is(x, "json")

    ess <- protolite::unserialize_pb(readBin(fromJSON(x), what="raw", n=100000000))

    expect_equal(length(ess), 1)

    expect_is(loadGEO("GSE14308"), "json")
    expect_is(loadGEO("GDS4885"), "json")

    expect_error(loadGEO("WRONGNAME"))

    options(phantasusMirrorPath = NULL)
})

test_that("reparseCachedGSEs works", {
    cacheDir <- tempdir()
    getES("GSE14308", destdir = cacheDir)
    expect_true("GSE14308" %in% reparseCachedESs(destdir = cacheDir))
})

test_that("checkGPLs counts gpls correctly", {
    options(phantasusMirrorPath = "https://genome.ifmo.ru/files/software/phantasus")

    expect_equal(fromJSON(checkGPLs("GSE14308")), c("GSE14308"))
    expect_equal(fromJSON(checkGPLs("GDS4885")), c("GDS4885"))
    expect_length(fromJSON(checkGPLs("GSE27112")), 2)
    expect_length(fromJSON(checkGPLs("GSE10000")), 2)
    expect_warning(checkGPLs("GSE101"))
    expect_warning(checkGPLs("GSE201"))

    options(phantasusMirrorPath = NULL)
})

test_that("checkGPLs counts existing files correctly without connection", {
    options(phantasusMirrorPath = "https://notworkingdomain",
            phantasusCacheDir = system.file("testdata", package="phantasus"))

    expect_message(checkGPLs("GSE27112"), regexp = "Problems establishing connection")
    expect_length(fromJSON(checkGPLs("GSE27112")), 1)
    expect_warning(checkGPLs("GSE14308"))

    options(phantasusCacheDir = NULL,
            phantasusMirrorPath = NULL)
})

test_that("getGSE works with ARCHS4", {
    ess <- getGSE("GSE99709", destdir=system.file("testdata", package="phantasus"))
    expect_gt(nrow(ess[[1]]), 0)
    expect_gt(ncol(ess[[1]]), 0)
})
