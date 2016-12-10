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
  ind <- which(is.na(data), arr.ind = T)
  if (nrow(ind) > 0) {
    data[ind] <- apply(data, 1, replacena, na.rm = T)[ind[,1]]
  }
  ind1 <- which(!is.nan(as.matrix(data)), arr.ind = T)
  left.rows <- unique(ind1[,"row"])
  data <- data[left.rows,]
  data <- t(scale(t(data)))
  km <- stats::kmeans(data, k)
  res <- data.frame(row.names = row.names(exprs(es)))
  res[["cluster"]] <- NA
  res[names(km$cluster), "cluster"] <- as.vector(km$cluster)
  return(toJSON(as.vector(res[["cluster"]])))
}

