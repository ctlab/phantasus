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
    data <- data[rows, ]
    data <- replacenas(data, replacena)
    data <- t(scale(t(data)))
  }
  data
}

replacenas <- function(data, replacena) {
    ind <- which(is.na(data), arr.ind = TRUE)
    if (nrow(ind) > 0) {
        data[ind] <- apply(data, 1, replacena, na.rm = TRUE)[ind[, 1]]
    }
    ind1 <- which(!is.nan(as.matrix(data)), arr.ind = TRUE)
    left.rows <- unique(ind1[, "row"])
    data <- data[left.rows, ]
    data
}

filternaRows <- function(data, currentRows) {
  sums <- rowSums(data)
  rows <- currentRows[!(currentRows %in% which(is.na(sums)))]
  rows
}

#' Reads ExpressionSet from GCT file.
#'
#' Only versions 1.2 and 1.3 are supported.
#'
#' @param gct Path to gct file
#'
#' @param ... additional options for read.csv
#'
#' @return ExpressionSet object
#' @export
#'
#' @examples
#' read.gct(system.file("extdata", "centers.gct", package = "phantasus"))
read.gct <- function(gct, ...) {
  meta <- readLines(gct, n = 3)
  version <- meta[1]
  size <- as.numeric(unlist(strsplit(meta[2], "\t")))

  if (grepl("^#1.3", version)) {
    # number of column annotations = number of additional rows
    ann.col <- size[4]

    # number of row annotations = number of additional columns
    ann.row <- size[3]
  } else if (grepl("^#1.2", version)) {
    ann.col <- 0
    ann.row <- 1
  } else {
    stop("Unsupported version of gct: use 1.2 or 1.3")
  }

  t <- read.tsv(gct, skip = 2 + 1 + ann.col, nrows = size[1],
                col.names = unlist(strsplit(meta[3], "\t")),
                row.names = 1, header = FALSE,  ...)

  exp <- as.matrix(t[, (ann.row + 1):ncol(t)])

  fdata <- makeAnnotated(t[, seq_len(ann.row), drop = FALSE])


  if (ann.col > 0) {
    pdata.raw <- t(read.tsv(gct, skip = 2 + 1, nrows = ann.col, header = FALSE))
    pdata <- data.frame(pdata.raw[seq_len(ncol(exp)) + 1 + ann.row, ,
                                  drop = FALSE])
    colnames(pdata) <- pdata.raw[1, ]
    rownames(pdata) <- colnames(exp)
    pdata <- makeAnnotated(pdata)

    res <- ExpressionSet(exp, featureData = fdata, phenoData = pdata)
  } else {
    res <- ExpressionSet(exp, featureData = fdata)
  }

  res
}

read.tsv <- function(file, header = TRUE, sep = "\t", quote = "",
                     comment.char = "",
                     check.names = FALSE, ...) {
  args <- list(...)
  res <- utils::read.table(file, header = header, sep = sep, quote = quote,
                    comment.char = comment.char, check.names = check.names,
                    stringsAsFactors = FALSE,
                    ...)
  if ( (!"row.names" %in% names(args)) && (colnames(res)[1] == "") ) {
    rownames(res) <- res[, 1]
    res[[1]] <- NULL
  }
  res
}

makeAnnotated <- function(data) {
  meta <- data.frame(labelDescription = colnames(data))
  rownames(meta) <- colnames(data)

  methods::new("AnnotatedDataFrame", data = data, varMeta = meta)
}
