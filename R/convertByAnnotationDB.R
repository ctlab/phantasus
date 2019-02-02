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
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#' queryAnnotationDBMeta()
#' }
#'
queryAnnotationDBMeta <- function () {
    cacheDir <- getOption("phantasusCacheDir")
    annotDir <- paste(cacheDir, "annotationdb", sep = .Platform$file.sep)
    columnFiles <- list.files(annotDir, '.selected_fields.txt', full.names = TRUE)
    metaList <- list()
    for(columnFile in columnFiles) {
        columnsTSV <- read.table(file = columnFile, sep = '\t', header = TRUE, skip = 1)
        columns <- apply(columnsTSV, 1, paste, collapse = " - ")

        orgName = unlist(strsplit(basename(columnFile), ".selected_fields.txt"))[1]

        metaList[[orgName]] = list(species = readLines(columnFile, n = 1), columns = columns)
    }

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

    message('Populating AnnotationDB cache')
    dbFiles <- list.files(annotDir, '*.sqlite$', full.names = TRUE)
    for (dbFile in dbFiles) {
        db <- loadDb(dbFile)
        columnFile <- paste(dbFile, ".selected_fields.txt", sep = "")
        if (!file.exists(columnFile)) {
            columnsDB <- columns(db)
            columns <- sapply(columnsDB, function (column) {
                hint <- paste(head(keys(db, keytype=column), n = 2), collapse=";")
                return(paste(column, hint, sep=" - "))
            })

            humanMeta <- t(data.frame(strsplit(columns, " - ", fixed=TRUE)))
            cat(species(db), "\n", file=columnFile)
            suppressWarnings(write.table(humanMeta, file=columnFile, quote=FALSE, sep='\t',  row.names=FALSE, col.names=c('FIELD', 'HINT'), append=TRUE))
        }
    }
}
