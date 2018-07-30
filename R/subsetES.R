#' Subsets es, if rows or columns are not specified, all are retained
#' @param es ExpressionSet object.#'
#' @param columns List of specified columns' indices (optional), indices start from 0#'
#' @param rows List of specified rows' indices (optional), indices start from 0
#' @return new expression set `es`
#'
#'
subsetES <- function(es, columns = c(), rows=c()) {
    rows <- getIndicesVector(rows, nrow(exprs(es)))
    columns <- getIndicesVector(columns, ncol(exprs(es)))

    es <- es[rows, columns]
    assign("es", es, envir = parent.frame())
}
