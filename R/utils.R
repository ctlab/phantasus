getIndicesVector <- function(current, neededLength) {
    if (length(current) == 0) {
        current <- 0:(neededLength - 1)
    }
    current + 1
}

prepareData <- function(es, columns = c(), rows = c(), replacena = "mean") {
  rows <- getIndicesVector(rows, nrow(exprs(es)))
  columns <- getIndicesVector(columns, ncol(exprs(es)))
  data <- replacenas(data.frame(exprs(es))[rows, columns], replacena)

  data <- t(scale(t(data)))
  while (sum(is.na(data)) > 0) {
    rows <- filternaRows(data, rows)
    data <- data[rows,]
    data <- replacenas(data, replacena)
    data <- t(scale(t(data)))
  }
  data
}

replacenas <- function(data, replacena) {
    ind <- which(is.na(data), arr.ind = T)
    if (nrow(ind) > 0) {
        data[ind] <- apply(data, 1, replacena, na.rm = T)[ind[, 1]]
    }
    ind1 <- which(!is.nan(as.matrix(data)), arr.ind = T)
    left.rows <- unique(ind1[, "row"])
    data <- data[left.rows, ]
    data
}

filternaRows <- function(data, currentRows) {
  sums <- rowSums(data)
  rows <- currentRows[!(currentRows %in% which(is.na(sums)))]
  rows
}
