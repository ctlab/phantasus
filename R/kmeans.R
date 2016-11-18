kmeans <- function(df, cols = с(), rows = с(), k, max.iterations = 10000) {
  if (!is.null(cols)) {
    df <- df[,cols]
  }
  if (!is.null(rows)) {
    df <- df[rows,]
  }
  used.df <- normalizeDF(df)
  used.df[["label"]] <- NA
  centroids <- used.df[sample(1:nrow(used.df), size = k),]
  centroids[["label"]] <- 1:k
  iterations <- 0
  old.centroids <- NULL
  while(!enough(old.centroids, centroids, iterations, max.iterations)) {
    old.centroids <- centroids
    iterations <- iterations + 1
    used.df <- labelDF(used.df[,1:(ncol(used.df) - 1)], centroids)
    centroids <- updateCentroids(used.df, k)
  }
  return(used.df[["label"]])
}

normalizeDF <- function(df) {
  res <- df
  for(i in 1:nrow(df)) {
    res[i,] <- (df[i,] - mean(as.matrix(df[i,]))) / sd(as.matrix(df[i,]))
  }
  return(res)
}

updateCentroids <- function(df, k) {
  centroids <- data.frame(matrix(0, k, ncol(df)))
  colnames(centroids) <- colnames(df)
  for(i in 1:k) {
    cur <- df[df$label == i, 1:(ncol(df) - 1)]
    centroids[i,] <- c(as.vector(colSums(cur) / nrow(cur)), i)
  }

  return(centroids)
}

labelDF <- function(df, centroids) {
  labels <- c()
  for(i in 1:nrow(df)) {
    labels <- c(labels, chooseClosestCentroid(df[i,], centroids))
  }
  df[["label"]] <- labels
  return(df)
}

enough <- function(old.centroids, centroids, iterations, max.iterations) {
  return(iterations > max.iterations || !is.null(old.centroids) && old.centroids == centroids)
}

chooseClosestCentroid <- function(x, centroids) {
  min <- 0
  res <- 0
  for(i in 1:nrow(centroids)) {
    d <- distance(x, centroids[i, 1:(ncol(centroids) - 1)])
    if (abs(d) > min) {
      res <- centroids[i, "label"]
      min <- d
    }
  }
  return(res)
}

distance <- function(x, y) {
  return(cor(c(as.matrix(x)), c(as.matrix(y)), method = "pearson"))
}

