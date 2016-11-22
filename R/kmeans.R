kmeans <- function(es, cols = c(), rows = c(), k) {
  stopifnot(require(Biobase))
  used.df <- data.frame(exprs(es))
  if (!is.null(cols)) {
    used.df <- used.df[,cols]
  }
  if (!is.null(rows)) {
    used.df <- used.df[rows,]
  }
  used.df <- t(scale(t(used.df)))
  km <- stats::kmeans(used.df, k)
  res <- data.frame(row.names = row.names(exprs(ALL)))
  res[["cluster"]] <- NA
  res[names(km$cluster), "cluster"] <- as.vector(km$cluster)
  return(as.vector(res[["cluster"]]))
}

