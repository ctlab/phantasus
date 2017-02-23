quantileNormalization <- function(es, cols = c(), rows = c()) {
  stopifnot(require(limma))
  stopifnot(require(Biobase))
  stopifnot(require(jsonlite))
  exprs <- exprs(es)
  if (is.null(cols)) {
    cols <- 0:(ncol(exprs) - 1)
  }
  if (is.null(rows)) {
    rows <- 0:(nrow(exprs) - 1)
  }
  exprs <- exprs[rows + 1, cols + 1]
  exprs(es)[rows, cols] <- normalizeBetweenArrays(log2(exprs + 1), method = "quantile")
  assign("es", es, envir = parent.frame())
  jsonlite::toJSON(exprs(es)[rows, cols])
}
