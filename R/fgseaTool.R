availableFGSEADatabases <- function () {
    dbDir <- getPhantasusConf("local_cache")$fgsea_pathways
    metaFile <- file.path(dbDir, 'fgsea.txt')
    meta = read.table(metaFile, header = TRUE, sep = '\t')

    return (jsonlite::toJSON(meta))
}

FGSEAmeta <- function (dbDir) {
    if (!dir.exists(dbDir)) {
        return()
    }

    metaFile <- file.path(dbDir, 'fgsea.txt')
    columnFiles <- list.files(dbDir, '*.rds', full.names = FALSE)
    message('Populating FGSEA pathway meta')
    if (file.exists(metaFile)) {
        meta = read.table(metaFile, header = TRUE, sep = '\t')
    } else {
        meta = data.frame(matrix(ncol=2,nrow=0), stringsAsFactors = FALSE)
    }
    colnames(meta) <- c('FILE', 'HINT')

    for(file in columnFiles) {
        if (nrow(meta[meta$FILE == file,]) == 0) {
            meta[nrow(meta) + 1, ] <- c(file, file)
        }
    }

    write.table(meta, file=metaFile, quote=FALSE, sep='\t', row.names = FALSE)
}

performFGSEA <- function (dbName, ranks) {
    dbDir <- getPhantasusConf("local_cache")$fgsea_pathways
    dbFile <- file.path(dbDir, dbName)
    if (!file.exists(dbFile)) {
        stop('Invalid DB name supplied')
    }

    db <- readRDS(dbFile)
    pathways <- lapply(lapply(split(db, db$pathName), '[[', 'geneID'), as.character)
    rranks <- ranks$ranks
    names(rranks) <- ranks$genes

    assign('db', db, envir = parent.frame())

    return (jsonlite::toJSON(fgseaMultilevel(pathways,
                                             rranks,
                                             sampleSize = 101,
                                             eps = 1e-10,
                                             minSize = 15,
                                             maxSize = 500,
                                             nproc = 1)
                             , digits = NA
                             , na="string"))
}

queryPathway <- function (dbName, pathwayName) {
    dbDir <- getPhantasusConf("local_cache")$fgsea_pathways
    dbFile <- file.path(dbDir, dbName)
    if (!file.exists(dbFile)) {
        stop('Invalid DB name supplied')
    }

    db <- readRDS(dbFile)
    filteredDb <- db[db$pathName == pathwayName,]

    jsonlite::toJSON(list(geneID = filteredDb$geneID, geneSymbol = filteredDb$geneSymbol))
}
