quantileNormalization <- function(es, cols = c(), rows = c()) {
    exprs <- exprs(es)
    if (is.null(cols)) {
        cols <- 0:(ncol(exprs) - 1)
    }
    if (is.null(rows)) {
        rows <- 0:(nrow(exprs) - 1)
    }
    exprs <- exprs[rows + 1, cols + 1]
    exprs(es)[rows + 1, cols + 1] <- normalizeBetweenArrays(log2(exprs + 1),
                                                            method = "quantile")
    assign("es", es, envir = parent.frame())
    f <- tempfile(pattern = "qn", tmpdir = getwd(), fileext = ".bin")
    writeBin(protolite::serialize_pb(list(data = exprs(es)[rows + 1, cols + 1])), f)
    f
}
