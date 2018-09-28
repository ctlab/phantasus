#' Create calculated annotation
#'
#' \code{calculatedAnnotation} adds a column calculated by operation
#'
#' @param es ExpressionSet object.
#' @param operation Name of the operation to perform calculation
#' @param columns List of specified columns' indices (optional), indices start from 0#'
#' @param rows List of specified rows' indices (optional), indices start from 0
#' @param name Name of the new annotation
#' @param isColumns Apply fn to columns
#'
#' @import Biobase
#'

calculatedAnnotation <- function (es, operation, rows = c(), columns = c(), isColumns = FALSE, name = NULL) {
    rows <- getIndicesVector(rows, nrow(exprs(es)))
    columns <- getIndicesVector(columns, ncol(exprs(es)))

    fn <- tolower(operation)
    if (!isColumns) {
        fData(es)[[name]] <- NA
        fData(es)[[name]][rows] <- apply(exprs(es[rows,columns]), 1, fn)
    } else {
        phenoData(es)[[name]] <- NA
        phenoData(es)[[name]][columns] <- apply(exprs(es[rows,columns]), 2, fn)
    }

    assign("es", es, envir = parent.frame())

}
