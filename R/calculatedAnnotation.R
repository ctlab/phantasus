#' Create calculated annotation
#'
#' \code{calculatedAnnotation} adds a column calculated by operation
#'
#' @param es ExpressionSet object.
#'
#' @param operation Name of the operation to perform calculation
#'
#' @import Biobase
#'

calculatedAnnotation <- function (es, operation) {
    fn <- tolower(operation)
    fData(es)[[operation]] <- apply(exprs(es), 1, fn)
    assign("es", es, envir = parent.frame())

}
