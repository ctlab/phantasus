context("Load GEO and its utils")
library(jsonlite)
library(Biobase)
library(data.table)

test_that("loadGEO finishes with result", {
    options(phantasusMirrorPath = "https://genome.ifmo.ru/files/software/phantasus",
            phantasusCacheDir = tempdir())

    cacheDir <- getOption("phantasusCacheDir")
    x <- loadGEO("GSE27112")
    expect_is(x, "json")

    binPath <- file.path(cacheDir, fromJSON(x))
    ess <- protolite::unserialize_pb(readBin(binPath, what="raw", n=100000000))

    expect_equal(length(ess), 2)

    x <- loadGEO("GSE27112-GPL6885")
    expect_is(x, "json")

    binPath <- file.path(cacheDir, fromJSON(x))
    ess <- protolite::unserialize_pb(readBin(binPath, what="raw", n=100000000))

    expect_equal(length(ess), 1)

    expect_is(loadGEO("GSE14308"), "json")
    expect_is(loadGEO("GDS4885"), "json")

    expect_error(loadGEO("WRONGNAME"))

    options(phantasusMirrorPath = NULL, phantasusCacheDir = NULL)
})

test_that("getGDS adds id field for GDS datasets", {
    a <- getGDS("GDS4885")[[1]]
    expect_true("id" %in% tolower(fvarLabels(a)))
})

test_that("filterPhenoAnnotations saves colnames", {
    cacheDir <- tempdir()
    es <- getES("GSE53986", destdir = cacheDir)[[1]]
    expect_true(all(colnames(es) == colnames(exprs(es))))
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

test_that("checkGPLs works with fully specified name", {
    options(phantasusMirrorPath = "https://genome.ifmo.ru/files/software/phantasus")

    expect_equal(fromJSON(checkGPLs("GSE27112-GPL6885")), c("GSE27112-GPL6885"))

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

test_that("InferConditionImpl  works correctly", {
    tests <- fread(system.file("testdata/dts.tsv", package="phantasus"))
    test_ds <- data.table(title=tests$Title, series=tests$Series, accession=tests$Accession, rep=tests$Replicate, inferCondition=tests$InferCondition)
    cond <- split(test_ds$title, test_ds$series)
    inf_cond_test <- split(test_ds$inferCondition, test_ds$series)
    rep_test <- split(test_ds$rep, test_ds$series)
    new_cond <- lapply(cond, inferConditionImpl)
    expect_equal(new_cond$GSE100221, list()) # text in all titles is unique
    expect_equal(new_cond$GSE10380, list())  # long dataset
    expect_equal(new_cond$GSE10382, list())  # number-only titles
    expect_equal(new_cond$GSE10383, list())  # two-color datasets
    expect_equal(new_cond$GSE10385, list())  # the same text and replicate number in all titles
    expect_equal(new_cond$GSE10039, list())  # ambiguous replicate number "High_Mo_seg_pool_Ler_col_F2" "Low_Mo_seg_pool_Ler_col_F2"  "Col-0 3"
    expect_equal(new_cond$GSE101508$condition, inf_cond_test$GSE101508) #"IFNγ+LPS rep2" -> "IFNγ+LPS" + "2"
    expect_equal(new_cond$GSE101508$replicate, as.character(rep_test$GSE101508))
    expect_equal(new_cond$GSE10392$condition, inf_cond_test$GSE10392) # "MPA 1" - > "MPA" + "1"
    expect_equal(new_cond$GSE10392$replicate, as.character(rep_test$GSE10392))
    expect_equal(new_cond$GSE10123$condition, inf_cond_test$GSE10123) # "WT-GFP-lamin A Induction: Day 0 Replicate A" -> "WT-GFP-lamin A Induction: Day 0" + "A"
    expect_equal(new_cond$GSE10123$replicate, as.character(rep_test$GSE10123))


})

