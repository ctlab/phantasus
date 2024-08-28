#' K-means clusterisation.
#'
#' \code{performKmeans} returns a vector of corresponding clusters for
#'     each gene from a given ExpressionSet.
#'
#' @param es ExpressionSet object.
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
#' @keywords internal
performKmeans <- function(es, k,replacena = "mean") {
    assertthat::assert_that(k > 0)

    scaledExprs <- unname(exprs(es))

    naInd <- which(is.na(scaledExprs), arr.ind = TRUE)
    if (nrow(naInd) > 0) {
        replaceValues <- apply(scaledExprs, 1, replacena, na.rm=TRUE)
        scaledExprs[naInd] <- replaceValues[naInd[,1]]
    }

    scaledExprs <- t(scale(t(scaledExprs)))
    rowsToCluster <- which(!apply(is.na(scaledExprs), 1, any))

    km <- stats::kmeans(scaledExprs[rowsToCluster, ], k, iter.max = 100L)
    res <- character(nrow(scaledExprs))
    res[rowsToCluster] <- as.character(km$cluster)
    fData(es)$clusters = res
    assign("es", es, envir = parent.frame())

    return(jsonlite::toJSON(res))
}

