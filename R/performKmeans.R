#' K-means clusterisation.
#'
#' \code{performKmeans} returns a vector of corresponding clusters for
#'     each gene from a given ExpressionSet.
#'
#' @param es ExpressionSet object.
#'
#' @param columns List of specified columns' indices (optional), indices start from 0
#'
#' @param rows List of specified rows' indices (optional), indices start from 0
#'
#' @param k Expected number of clusters.
#'
#' @param replacena Method for replacing NA values
#'     in series matrix (mean by default)
#'
#' @return Vector of corresponding clusters, serialized to JSON.
#'
#' @import Biobase
#'
#' @examples
#' \dontrun{
#' data(es)
#' performKmeans(es, k = 2)
#' }
performKmeans <- function(es, columns = c(), rows = c(), k,
                            replacena = "mean") {
    assertthat::assert_that(k > 0)

    es <- subsetES(es, columns=columns, rows=rows)

    scaledExprs <- unname(t(scale(t(exprs(es)))))



    naInd <- which(is.na(scaledExprs), arr.ind = TRUE)
    if (nrow(naInd) > 0) {
        replaceValues <- apply(scaledExprs, 1, replacena, na.rm=T)
        scaledExprs[naInd] <- replaceValues[naInd[,1]]
        rowsToCluster <- is.finite(replaceValues)
    } else {
        rowsToCluster <- seq_len(nrow(scaledExprs))
    }

    km <- stats::kmeans(scaledExprs[rowsToCluster, ], k, iter.max = 100L)
    res <- character(nrow(scaledExprs))
    res[rowsToCluster] <- as.character(km$cluster)
    return(jsonlite::toJSON(res))
}
