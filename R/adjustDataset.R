adjustDataset <- function (es, scaleColumnSum = NULL, log2 = FALSE, onePlusLog2 = FALSE,
                           inverseLog2 = FALSE, quantileNormalize = FALSE, zScore = FALSE,
                           robustZScore = FALSE, sweep = NULL) {

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

