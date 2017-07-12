#' Principal Component Analysis.
#'
#' Function for creating json with full description
#'   of the pcaPlot for plotly.js.
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
#' @export
#' @import ggplot2
#' @import htmltools
#' @import jsonlite
pcaPlot <- function(es, columns = c(), rows = c(), replacena = "mean") {

    data <- prepareData(es, columns, rows, replacena)

    pca <- stats::prcomp(data)
    explained <- (pca$sdev) ^ 2 / sum(pca$sdev ^ 2)

    xs <- sprintf("PC%s", seq_along(explained))
    xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)

    pca.res <- as.matrix(pca$x)
    colnames(pca.res) <- NULL
    row.names(pca.res) <- NULL
    return(jsonlite::toJSON(list(pca = t(pca.res), xlabs = xlabs)))
}
