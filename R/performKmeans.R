#' K-means clusterisation.
#'
#' \code{performKmeans} returns a vector of corresponding clusters for
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
#'
#' @examples
#' data(es)
#' performKmeans(es, k = 2)
#'
performKmeans <- function(es, columns = c(), rows = c(), k,
                          replacena = "mean") {
    assertthat::assert_that(k > 0)

    data <- prepareData(es, columns, rows, replacena)

    km <- stats::kmeans(data, k, iter.max = 100L)
    res <- data.frame(row.names = row.names(exprs(es)))
    res[["cluster"]] <- NA
    res[names(km$cluster), "cluster"] <- as.vector(km$cluster)
    return(jsonlite::toJSON(as.vector(km$cluster)))
}
