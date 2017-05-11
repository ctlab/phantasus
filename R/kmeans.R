#' @name kmeans
#' @title KMeans
#' @description Function for performing k-means clusterisation on ExpressionSet.
#' @param es an ExpressionSet object
#' @param columns list of specified columns' indices (optional)
#' @param rows list of specified rows' indices (optional)
#' @param k expected number of clusters
#' @param replacena method for replacing NA values (mean by default)
#' @return json, containing array of corresponding clusters
#' @export
#' @import Biobase
#' @import jsonlite
#' @import assertthat
kmeans <- function(es, columns = c(), rows = c(), k, replacena = "mean") {
  assertthat::assert_that(k > 0)

  rows <- getIndicesVector(rows, nrow(exprs(es)))
  columns <- getIndicesVector(columns, ncol(exprs(es)))
  data <- replacenas(data.frame(exprs(es))[rows, columns], replacena)

  data <- t(scale(t(data)))
  while (sum(is.na(data)) > 0) {
    data <- replacenas()
    data <- t(scale(t(data)))
  }

  km <- stats::kmeans(data, k)
  res <- data.frame(row.names = row.names(exprs(es)))
  res[["cluster"]] <- NA
  res[names(km$cluster), "cluster"] <- as.vector(km$cluster)
  return(toJSON(as.vector(km$cluster)))
}

