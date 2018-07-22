getIndicesVector <- function(current, neededLength) {
    if (length(current) == 0) {
        current <- 0:(neededLength - 1)
    }
    current + 1
}

#' Subsets es, if rows or columns are not specified, all are retained
#' @param es ExpressionSet object.#'
#' @param columns List of specified columns' indices (optional), indices start from 0#'
#' @param rows List of specified rows' indices (optional), indices start from 0
#' @return `es`'s subset
#'
subsetES <- function(es, columns = c(), rows=c()) {
    rows <- getIndicesVector(rows, nrow(exprs(es)))
    columns <- getIndicesVector(columns, ncol(exprs(es)))

    es[rows, columns]
}

prepareData <- function(es, columns = c(), rows = c(), replacena = "mean") {
    rows <- getIndicesVector(rows, nrow(exprs(es)))
    columns <- getIndicesVector(columns, ncol(exprs(es)))

    data <- replacenas(data.frame(exprs(es[rows, columns])), replacena)

    rows <- getIndicesVector(c(), nrow(data))

    data <- t(scale(t(data)))
    while (sum(is.na(data)) > 0) {
        message("need to filter rows")
        rows <- filternaRows(data, rows)

        message(length(rows))
        message(rows[length(rows)])

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

#' Reads ExpressionSet from a GCT file.
#'
#' Only versions 1.2 and 1.3 are supported.
#'
#' @param gct Path to gct file
#'
#' @param ... additional options for read.csv
#'
#' @return ExpressionSet object
#'
#' @examples
#' read.gct(system.file("extdata", "centers.gct", package = "phantasus"))
#' @export
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

    colNames <- unlist(strsplit(meta[3], "\t"))
    if (grepl("/", colNames[1])) {
        rowIdField <- sub("(.*)/(.*)", "\\1", colNames[1])
        colIdField <- sub("(.*)/(.*)", "\\2", colNames[1])
    } else {
        rowIdField <- "id"
        colIdField <- "id"
    }

    colNames[1] <- rowIdField

    t <- read.tsv(gct, skip = 2 + 1 + ann.col, nrows = size[1],
                  col.names = colNames,
                  row.names = NULL, header = FALSE,  ...)

    if (any(duplicated(t[,1]))) {
        warning(sprintf("duplicated row IDs: %s; they were renamed",
                        paste0(t[head(which(duplicated(t[, 1]))),1], collapse = " ")))
        rownames(t) <- make.unique(t[,1])
    } else {
        rownames(t) <- t[,1]
    }


    exp <- as.matrix(t[, (ann.row + 2):ncol(t)])

    fdata <- makeAnnotated(t[, seq_len(ann.row + 1), drop = FALSE])


    if (ann.col > 0) {
        pdata.raw <- t(read.tsv(gct, skip = 2, nrows = ann.col + 1,
                                header = FALSE, row.names=NULL))
        pdata <- data.frame(pdata.raw[seq_len(ncol(exp)) + 1 + ann.row, ,
                                      drop = FALSE])
        colnames(pdata) <- pdata.raw[1, ]
        colnames(pdata)[1] <- colIdField
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

#' Saves ExpressionSet to a GCT file (version 1.3).
#'
#' @param es ExpresionSet obeject to save
#' @param file Path to output gct file
#' @param gzip Whether to gzip apply gzip-compression for the output file#'
#' @param ... additional options for read.csv
#' @return Result of the closing file (as in `close()` function`)
#' @examples
#' es <- read.gct(system.file("extdata", "centers.gct", package = "phantasus"))
#' out <- tempfile(fileext = ".gct.gz")
#' write.gct(es, out, gzip=TRUE)
#' @import Biobase
#' @export
write.gct <- function(es, file, gzip=FALSE) {
    if (gzip) {
        con <- gzfile(file)
    } else {
        con <- file(file)
    }
    open(con, open="w")
    writeLines("#1.3", con)
    ann.col <- ncol(pData(es))
    ann.row <- ncol(fData(es))
    writeLines(sprintf("%s\t%s\t%s\t%s", nrow(es), ncol(es), ann.row, ann.col), con)
    writeLines(paste0(c("ID", colnames(fData(es)), colnames(es)), collapse="\t"), con)

    ann.col.table <- t(as.matrix(pData(es)))
    ann.col.table <- cbind(matrix(rep(NA, ann.row*ann.col), nrow=ann.col), ann.col.table)
    write.table(ann.col.table, file=con, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
    write.table(cbind(fData(es), exprs(es)), file=con, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
    close(con)
}


makeAnnotated <- function(data) {
    meta <- data.frame(labelDescription = colnames(data))
    rownames(meta) <- colnames(data)

    methods::new("AnnotatedDataFrame", data = data, varMeta = meta)
}

take <- function(x, n) {
  sapply(x, function(x) {
    x[[n]]
  })
}

writeToList <- function(es) {
  data <- as.matrix(exprs(es))
  colnames(data) <- NULL
  row.names(data) <- NULL

  pdata <- as.matrix(pData(es))
  colnames(pdata) <- NULL
  row.names(pdata) <- NULL

  rownames <- rownames(es)

  fdata <- as.matrix(fData(es))
  colnames(fdata) <- NULL
  row.names(fdata) <- NULL

  res <- list(data = data, pdata = pdata, fdata = fdata,
              rownames = rownames,
              colMetaNames = varLabels(es),
              rowMetaNames = fvarLabels(es))
  res
}
