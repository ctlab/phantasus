#' Adjust dataset
#'
#' @param es Expression set to perform adjustment on
#'
#' @param scaleColumnSum perform sum scaling of columns (default FALSE)
#'
#' @param log2 perform logarithm2 adjustment (default FALSE)
#'
#' @param onePlusLog2 perform log2(1+x) adjustment (default FALSE)
#'
#' @param inverseLog2 perform 2^x adjustment (default FALSE)
#'
#' @param quantileNormalize perform quantile normalization (default FALSE)
#'
#' @param zScore perform zScore adjustment:
#'  subtract mean, divide by std (default FALSE)
#'
#' @param robustZScore perform robustZScore adjustment:
#'  subtract median, divide by MAD (default FALSE)
#'
#' @param sweep perform sweep adjustment on rows/columns (default FALSE)
#'
#' @return Nothing. Adjusted dataset will be assigned as ES in global environment
#'
#' @import stats
#'
#' @examples
#' \dontrun{
#' es <- gseGSE('GSE53986')[[1]]
#' adjustDataset(es, log2 = T, quantileNormalize = T)
#' }
#'

adjustDataset <- function (es, scaleColumnSum = NULL,
                           log2 = FALSE, onePlusLog2 = FALSE,
                           inverseLog2 = FALSE, quantileNormalize = FALSE,
                           zScore = FALSE, robustZScore = FALSE, sweep = NULL) {

    if (length(scaleColumnSum) > 0 && is.numeric(scaleColumnSum)) {
        sum <- apply(exprs(es), 2, sum, na.rm = TRUE)
        sum <- t(scaleColumnSum / sum)
        exprs(es) <- sweep(exprs(es), 2, sum, "*")
    }

    if (log2) {
        exprs(es) <- apply(exprs(es), c(1,2), function (value) {
            safeLog2(value)
        })
    }

    if (onePlusLog2) {
        exprs(es) <- apply(exprs(es), c(1,2), function (value) {
            safeLog2(value + 1)
        })
    }

    if (inverseLog2) {
        exprs(es) <- apply(exprs(es), c(1,2), function (value) {
            2^value
        })
    }

    if (quantileNormalize) {
        exprs(es) <- normalizeBetweenArrays(exprs(es), method = 'quantile')
    }

    if (zScore) {
        exprs(es) <- t(scale(t(exprs(es))))
    }

    if (robustZScore) {
        medians <- apply(exprs(es), 1, median)
        exprs(es) <- sweep(exprs(es), 1, medians, "-")
        mads <- apply(exprs(es), 1, mad)
        exprs(es) <- sweep(exprs(es), 1, mads, "/")
    }

    if (length(sweep) > 0) {
        if (sweep$mode == 'row') {
            margin <- 1
            target <- fData(es)[[sweep$name]]
        } else {
            margin <- 2
            target <- phenoData(es)[[sweep$name]]
        }

        exprs(es) <- sweep(exprs(es), margin, target, sweep$op)
    }

    assign("es", es, envir = parent.frame())
}

