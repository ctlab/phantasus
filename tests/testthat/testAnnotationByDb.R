context("Convert by AnnotationDb")
library(Biobase)
library(data.table)
test_that("AnnotationbyDb works with delete version", {
    options(phantasusCacheDir = system.file("testdata/", package = "phantasus"))
    test_file <- system.file("testdata/counts_versioned_ids.gct", package="phantasus")
    if (test_file == ""){
        stop("test counts file doesn't exists")
    }
    es <- read.gct(test_file)
    dbName <- "sample_mouse_db.sqlite"
    columnName <- "id"
    columnType <- "ENSEMBL"
    keyType <- "SYMBOL"
    otherOptions <- list(deleteDotVersion = FALSE)
    expect_error(convertByAnnotationDB(es, dbName, columnName, columnType, keyType, otherOptions))
    otherOptions <- list(deleteDotVersion = TRUE)
    symbols <- jsonlite::fromJSON(convertByAnnotationDB(es, dbName, columnName, columnType, keyType, otherOptions))
    expect_gt(sum(!is.na(symbols)),0)

})
