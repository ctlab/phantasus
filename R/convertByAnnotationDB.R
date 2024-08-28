#' Map indexes using Annotation DB
#'
#' \code{convertByAnnotationDB} function returns \code{keyType} ids from \code{dbName} mapped to \code{columnName} in \code{es}.
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
#' @param otherOptions additional parameters for conversion. Currently only named boolean deleteDotVersion is not ignored.
#'
#' @return JSON object with a vector of converted IDs
#' @importFrom AnnotationDbi loadDb mapIds
#'
#' @keywords internal
convertByAnnotationDB <- function (es, dbName, columnName, columnType, keyType, otherOptions) {
    annotDir <- getPhantasusConf("cache_folders")$annot_db
    dbPath <- file.path(annotDir, dbName)
    if (!file.exists(dbPath)) {
        stop('Invalid database specified')
    }

    dbFile <- loadDb(dbPath)

    inputData <- fData(es)[[columnName]]
    if ("deleteDotVersion" %in% names(otherOptions)){
        if(otherOptions$deleteDotVersion){
            inputData <- sub(pattern = "\\.[0-9]+$", replacement = "", x = inputData)
        }
    }
    convertedData <- mapIds(dbFile,
                           keys = inputData,
                           keytype = columnType,
                           column = keyType,
                           multiVals = function (x) { paste0(x, collapse="///")})

    convertedData[ convertedData == 'NA' ] <- NA # AnnotationDB can produce 'NA' as string and <NA> as NA. Confusing
    fData(es)[[keyType]] <- convertedData
    assign("es", es, envir = parent.frame())

    return(jsonlite::toJSON(convertedData))
}

#' Get meta list for annotationDB files
#'
#' \code{queryAnnotationDBMeta} Function reads txt meta files for provided
#' sqlite annotation databases.
#'
#' @return meta info in JSON
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#' queryAnnotationDBMeta()
#' }
#' @keywords internal
queryAnnotationDBMeta <- function () {
    annotDir <- getPhantasusConf("cache_folders")$annot_db
    if (is.null(annotDir)){
        stop("Current Phantasus configuration doesn't support annotation databases. 'annot_db' setting for cache folders is missed.")
    }
    if (!dir.exists(annotDir)){
        stop("Current Phantasus configuration doesn't support annotation databases. AnnotationDB folder does not exist.")
    }
    columnFiles <- list.files(annotDir, '.selected_fields.txt', full.names = TRUE)
    dbFiles <- list.files(annotDir, '*.sqlite$', full.names = TRUE)
    if (length(dbFiles) != length(columnFiles)){
        stop("AnnotationDB folder was not properly configured. Some meta files are missed.")
    }
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
#' \code{annotationDBMeta} function creates txt files containing meta information of
#' provided sqlite files for AnnotationDB.
#'
#' @param annotDir path to folder with annotationDB sqlite files
#'
#' @return nothing
#'
#' @importFrom AnnotationDbi keys species columns loadDb
#'
#' @examples
#' \dontrun{
#' annotationDBMeta('/var/phantasus/cache')
#' }
#' @keywords internal
annotationDBMeta <- function (annotDir) {
    if (!dir.exists(annotDir)) {
        return()
    }

    message('Populating AnnotationDB meta')
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
            cat(AnnotationDbi::species(db), "\n", file=columnFile)
            suppressWarnings(write.table(humanMeta, file=columnFile, quote=FALSE, sep='\t',  row.names=FALSE, col.names=c('FIELD', 'HINT'), append=TRUE))
        }
    }
}
