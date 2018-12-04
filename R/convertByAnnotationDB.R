#' Map indexes using Annotation DB
#'
#' \code{createES} function creates an rds file containing meta information of
#' provided sqlite files for AnnotationDB
#'
#' @param es source ExpressionSet
#'
#' @param dbName name of AnnotationDB file
#'
#' @param columnName name of column in featureData of source ExpressionSet
#'
#' @param columnType Type of indexes in columnName
#'
#' @param keyType Type of mapped indexes
#'
#' @import AnnotationDbi
#'
#' @examples
#' \dontrun{
#' }
#'
#'
convertByAnnotationDB <- function (es, dbName, columnName, columnType, keyType) {
    cacheDir <- getOption("phantasusCacheDir")

    annotDir <- paste(cacheDir, "annotationdb", sep = .Platform$file.sep)
    dbPath <- file.path(annotDir, dbName)
    if (!file.exists(dbPath)) {
        stop('Invalid database specified')
    }

    dbFile <- loadDb(dbPath)

    inputData <- fData(es)[[columnName]]
    convertedData <- mapIds(dbFile,
                           keys = inputData,
                           keytype = columnType,
                           column = keyType,
                           multiVals = function (x) { paste0(x, collapse="///")})

    if (keyType == 'ENSEMBL') {
        # remove versions
    }

    fData(es)[[keyType]] <- convertedData
    assign("es", es, envir = parent.frame())

    return(jsonlite::toJSON(convertedData))
}

#' Get meta list for annotationDB files
#'
#' \code{createES} Function reads an rds file containing meta information of provded
#' sqlite files for AnnotationDB
#'
#' @return meta info in JSON
#'
#' @examples
#' \dontrun{
#' queryAnnotationDBMeta
#' }
#'
queryAnnotationDBMeta <- function () {
    cacheDir <- getOption("phantasusCacheDir")
    annotDir <- paste(cacheDir, "annotationdb", sep = .Platform$file.sep)
    metaFile <- file.path(annotDir, "meta.rds")
    if (!file.exists(metaFile)) {
        return("[]")
    }

    metaList <- readRDS(metaFile)

    return(jsonlite::toJSON(metaList))
}

#' Create meta file for AnnotationDB
#'
#' \code{createES} function creates an rds file containing meta information of
#' provided sqlite files for AnnotationDB
#'
#' @param cacheDir cacheDir for phantasus
#'
#' @import AnnotationDbi
#'
#' @examples
#' \dontrun{
#' annotationDBMeta('/var/phantasus/cache')
#' }
#'
annotationDBMeta <- function (cacheDir) {
    annotDir <- file.path(cacheDir, "annotationdb")
    if (!dir.exists(annotDir)) {
        message('No annotationdb files provided')
        return()
    }

    metaList <- list()
    dbFiles <- list.files(annotDir, '*.sqlite$', full.names = TRUE)
    for (dbFile in dbFiles) {
        db <- loadDb(dbFile)
        columnFile <- paste(dbFile, ".selected_fields.txt", sep = "")
        if (file.exists(columnFile)) {
            columnsTSV <- read.table(file = columnFile, sep = '\t', header = TRUE)
            columns <- apply(columnsTSV, 1, paste, collapse = " - ")
            metaList[[basename(dbFile)]] <- list(species=species(db), columns=columns)
        } else {
            columnsDB <- columns(db)
            columns <- sapply(columnsDB, function (column) {
                hint <- paste(head(keys(db, keytype=column), n = 2), collapse=";")
                return(paste(column, hint, sep=" - "))
            })

            metaList[[basename(dbFile)]] <- list(species=species(db), columns=columns)
        }
    }

    saveRDS(metaList, file = file.path(annotDir, "meta.rds"))
}
