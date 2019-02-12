availableFGSEADatabases <- function () {
    cacheDir <- getOption("phantasusCacheDir")

    dbDir <- file.path(cacheDir, 'fgsea')
    metaFile <- file.path(dbDir, 'fgsea.txt')
    meta = read.table(metaFile, header = TRUE, sep = '\t')

    return (jsonlite::toJSON(meta))
}

FGSEAmeta <- function (cacheDir) {
    dbDir <- file.path(cacheDir, 'fgsea')
    if (!dir.exists(dbDir)) {
        message('No fgsea files provided')
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
    cacheDir <- getOption("phantasusCacheDir")

    dbDir <- file.path(cacheDir, 'fgsea')
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
                                             absEps = 0.5e-10,
                                             minSize = 15,
                                             maxSize = 500,
                                             nproc = 1)
                             , digits = NA))
}

queryPathway <- function (dbName, pathwayName) {
    cacheDir <- getOption("phantasusCacheDir")

    dbDir <- file.path(cacheDir, 'fgsea')
    dbFile <- file.path(dbDir, dbName)
    if (!file.exists(dbFile)) {
        stop('Invalid DB name supplied')
    }

    db <- readRDS(dbFile)
    filteredDb <- db[db$pathName == pathwayName,]

    jsonlite::toJSON(list(geneID = filteredDb$geneID, geneSymbol = filteredDb$geneSymbol))
}
