#' K-means clusterisation.
#'
#' \code{kmeans} returns a vector of corresponding clusters for
#'   each gene from a given ExpressionSet.
#'
#' @param es ExpressionSet object.
#'
#' @param columns List of specified columns' indices (optional).
#'
#' @param rows List of specified rows' indices (optional).
#'
#' @param k Expected number of clusters.
#'
#' @param replacena Method for replacing NA values
#'   in series matrix (mean by default)
#'
#' @return Vector of corresponding clusters, serialized to JSON.
#'
#' @export
#' @import Biobase
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

    km <- stats::kmeans(data, k, iter.max = 100L)
    res <- data.frame(row.names = row.names(exprs(es)))
    res[["cluster"]] <- NA
    res[names(km$cluster), "cluster"] <- as.vector(km$cluster)
    return(jsonlite::toJSON(as.vector(km$cluster)))
}
