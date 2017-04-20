#' @name pcaPlot
#' @title PCA Plot
#' @description Function for creating json with full description of the pcaPlot for plotly.js
#' @param es an ExpressionSet object, should be normalized
#' @param columns list of specified columns' indices (optional)
#' @param rows list of specified rows' indices (optional)
#' @param replacena method for replacing NA values (mean by default)
#' @return json with full description of the plot for plotly.js
#' @export
#' @import ggplot2
#' @import ggrepel
#' @import Biobase
#' @import svglite
#' @import plotly
#' @import htmltools
#' @import jsonlite
pcaPlot <- function(es, rows=c(), columns = c(), replacena = "mean") {
  rows <- getIndicesVector(rows, nrow(exprs(es)))
  columns <- getIndicesVector(columns, ncol(exprs(es)))
  data <- data.frame(exprs(es))[rows, columns]

  ind <- which(is.na(data), arr.ind = T)
  if (nrow(ind) > 0) {
    data[ind] <- apply(data, 1, replacena, na.rm = T)[ind[,1]]
  }
  ind1 <- which(!is.nan(as.matrix(data)), arr.ind = T)
  left.rows <- unique(ind1[,"row"])
  data <- data[left.rows,]
  data <- t(data)

  pca <- prcomp(data)
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)

  xs <- sprintf("PC%s", seq_along(explained))
  xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)

  pca.res <- as.matrix(pca$x); colnames(pca.res) <- NULL; row.names(pca.res) <- NULL
  return(jsonlite::toJSON(list(pca = t(pca.res), xlabs=xlabs)))
}

