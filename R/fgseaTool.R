availableFGSEADatabases <- function () {
    cacheDir <- getOption("phantasusCacheDir")

    dbDir <- file.path(cacheDir, 'fgsea')
    columnFiles <- list.files(dbDir, full.names = FALSE)
    columnFiles <- sapply(strsplit(as.character(columnFiles), "\\."), "[[", 1)

    return (jsonlite::toJSON(columnFiles))
}

performFGSEA <- function (dbName, ranks) {
    cacheDir <- getOption("phantasusCacheDir")

    dbDir <- file.path(cacheDir, 'fgsea')
    db <- readRDS(file.path(dbDir, paste(dbName, '.rds', sep='')))
    pathways <- lapply(lapply(split(db, db$pathName), '[[', 'geneID'), as.character)
    rranks <- ranks$ranks
    names(rranks) <- ranks$genes

    assign('db', db, envir = parent.frame())

    return (jsonlite::toJSON(fgseaMultilevel(pathways,
                                             rranks,
                                             sampleSize = 101,
                                             absEps = 0.5e-10,
                                             minSize = 15,
                                             maxSize = 500,
                                             nproc = 1)))
}

queryPathway <- function (dbName, pathwayName) {
    cacheDir <- getOption("phantasusCacheDir")

    dbDir <- file.path(cacheDir, 'fgsea')
    db <- readRDS(file.path(dbDir, paste(dbName, '.rds', sep='')))
    filteredDb <- db[db$pathName == pathwayName,]

    jsonlite::toJSON(list(geneID = filteredDb$geneID, geneSymbol = filteredDb$geneSymbol))
}
