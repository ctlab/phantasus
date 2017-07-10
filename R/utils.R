getIndicesVector <- function(current, neededLength) {
    if (length(current) == 0) {
        current <- 0:(neededLength - 1)
    }
    current + 1
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

