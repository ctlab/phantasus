kmeans <- function(es, cols = c(), rows = c(), k) {
  stopifnot(require(Biobase))
  stopifnot(require(jsonlite))
  used.df <- data.frame(exprs(es))
  if (!is.null(cols)) {
    used.df <- used.df[,(cols + 1)]
  }
  if (!is.null(rows)) {
    used.df <- used.df[(rows + 1),]
  }
  used.df <- t(scale(t(used.df)))
  km <- stats::kmeans(used.df, k)
  res <- data.frame(row.names = row.names(exprs(es)))
  res[["cluster"]] <- NA
  res[names(km$cluster), "cluster"] <- as.vector(km$cluster)
  return(toJSON(as.vector(res[["cluster"]])))
}

