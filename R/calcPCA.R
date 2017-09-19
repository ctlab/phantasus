#' Principal Component Analysis.
#'
#' \code{calcPCA} calculates PCA-matrix for the given ExpressionSet
#'     and returns this matrix encoded to JSON.
#'
#' @param es an ExpressionSet object, should be normalized
#'
#' @param columns list of specified columns' indices (optional)
#'
#' @param rows list of specified rows' indices (optional)
#'
#' @param replacena method for replacing NA values (mean by default)
#'
#' @return json with full description of the plot for plotly.js
#'
#' @import ggplot2
#' @import htmltools
#' @import jsonlite
#'
#' @examples
#' \dontrun{
#' data(es)
#' calcPCA(es)
#' }
calcPCA <- function(es, columns = c(), rows = c(), replacena = "mean") {

    data <- t(prepareData(es, columns, rows, replacena))

    pca <- stats::prcomp(data)
    explained <- (pca$sdev) ^ 2 / sum(pca$sdev ^ 2)

    xs <- sprintf("PC%s", seq_along(explained))
    xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)

    pca.res <- as.matrix(pca$x)
    colnames(pca.res) <- NULL
    row.names(pca.res) <- NULL
    return(jsonlite::toJSON(list(pca = t(pca.res), xlabs = xlabs)))
}
