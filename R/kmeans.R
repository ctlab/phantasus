kmeans <- function(es, cols = c(), rows = c(), k, replacena = "mean") {
  stopifnot(require(Biobase))
  stopifnot(require(jsonlite))
  data <- data.frame(exprs(es))
  if (!is.null(cols)) {
    data <- data[,(cols + 1)]
  }
  if (!is.null(rows)) {
    data <- data[(rows + 1),]
  }
  for(i in 1:nrow(data)) {
    data[i,] <- replace(data[i,], is.na(data[i,]), do.call(replacena, list(x = as.matrix(data[i,]), na.rm = TRUE)))
  }
  data <- t(scale(t(data)))
  km <- stats::kmeans(data, k)
  res <- data.frame(row.names = row.names(exprs(es)))
  res[["cluster"]] <- NA
  res[names(km$cluster), "cluster"] <- as.vector(km$cluster)
  return(toJSON(as.vector(res[["cluster"]])))
}

