getIndicesVector <- function(current, neededLength) {
    if (length(current) == 0) {
        current <- 0:(neededLength - 1)
    }
    current + 1
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

    rownames(t) <- t[,1]

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

  ed <- experimentData(es)
  experimentList <- as.list(expinfo(ed))
  experimentList$other <- as.list(ed@other)
  experimentList$pubMedIds <- pubMedIds(ed)

  res <- list(data = data, pdata = pdata, fdata = fdata,
              rownames = rownames,
              colMetaNames = varLabels(es),
              rowMetaNames = fvarLabels(es),
              experimentData = experimentList)
  res
}

#' @importFrom utils download.file
updateARCHS4 <- function (cacheDir = "/var/phantasus/cache/archs4") {
    download.file(url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5",
                  destfile = paste(cacheDir, "human_matrix.h5", sep=.Platform$file.sep),
                  mode = "wb")

    download.file(url = "https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5",
                  destfile = paste(cacheDir, "mouse_matrix.h5", sep=.Platform$file.sep),
                  mode = "wb")
}

selfCheck <- function (cacheDir=getOption("phantasusCacheDir"),
                       preloadedDir=getOption("phantasusPreloadedDir"),
                       verbose=FALSE) {

    if (!is.null(preloadedDir) && dir.exists(preloadedDir)) {
        preloadedFiles <- list.files(preloadedDir, pattern = "\\.(gct|rda)$")
        message(paste(length(preloadedFiles), 'preloaded datasets are available'))
        if (verbose) {
            message(paste0(preloadedFiles, collapse=" "))
            message(" ")
        }
    } else {
        message('!!! Preloaded dir is not set')
    }

    archs4Files <- list.files(file.path(cacheDir, "archs4"),
                              pattern='\\.h5$')
    if (length(archs4Files)) {
        message(paste(length(archs4Files), 'archs4 files are available'))
        if (verbose) {
            message(paste0(archs4Files, collapse=" "))
            message(" ")
        }
    } else {
        message('!!! No archs4 provided RNA-seq will load without matrices')
    }

    annotDir <- file.path(cacheDir, "annotationdb")
    dbFiles <- list.files(annotDir, pattern='\\.sqlite$')
    if (length(dbFiles)) {
        message(paste(length(dbFiles), 'annotationDb are available'))
        if (verbose) {
            message(paste0(dbFiles, collapse=" "))
            message(" ")
        }
    } else {
        message('!!! No annotationDb provided')
    }

    fgseaDir <- file.path(cacheDir, 'fgsea')
    fgseaFiles <- list.files(fgseaDir, '\\.rds$', full.names = FALSE)
    if (length(fgseaFiles)) {
        message(paste(length(fgseaFiles), 'fgsea tables are available'))
        if (verbose) {
            message(paste0(fgseaFiles, collapse=" "))
            message(" ")
        }
    } else {
        message('!!! No fgsea tables provided')
    }
}

safeDownload <- function (url, dir, filename) {
  dest <- file.path(dir, filename)
  if (file.exists(dest)) {
    return()
  }

  tempDest <- tempfile(paste0(filename, ".load"), tmpdir=dir)
  utils::download.file(url, destfile = tempDestFile)
  file.rename(tempDest, dest)
}


getGEODir <- function (name, destdir = tempdir()) {
  type <- substr(name, 1, 3)
  GEO <- unlist(strsplit(name, "-"))[1]
  stub <- gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  gdsDirPath <- "%s/geo/datasets/%s/%s/soft"
  gseDirPath <- "%s/geo/series/%s/%s/matrix"
  if (type == 'GSE') {
    fullGEODirPath <- file.path(sprintf(gseDirPath, destdir, stub, GEO))
  } else {
    fullGEODirPath <- file.path(sprintf(gdsDirPath, destdir, stub, GEO))
  }
  dir.create(fullGEODirPath, showWarnings = FALSE, recursive = TRUE)

  fullGEODirPath
}

getBriefData <- function (name, destdir = tempdir()) {
  GEO <- unlist(strsplit(name, "-"))[1]
  GEOdir <- dirname(getGEODir(GEO, destdir))
  briefFile <- file.path(GEOdir, 'brief')
  if (file.exists(briefFile)) {
    message('Using cached brief file: ', briefFile)
    return (parseBriefData(readLines(briefFile)))
  }

  url <- sprintf("www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief", GEO)
  resp <- httr::GET(url)
  text <- httr::content(resp, "text", "utf-8")
  check <- grep('Could not', text)
  if (length(check)) {
    message('No such dataset: ', name)
    unlink(GEOdir, recursive = TRUE, force = TRUE)
    stop('Failed to download brief data on: ', GEO, '. No such dataset')
  } else {
    writeLines(text, briefFile)
    message('Stored brief data of ', GEO, ' at ', briefFile)
  }

  return (parseBriefData(readLines(briefFile)))
}

parseBriefData <- function(txt) {
  tmp <- txt[grep("!\\w*?_",txt)]
  tmp <- gsub("!\\w*?_",'',tmp)
  first.eq <- regexpr(' = ',tmp)
  tmp <- cbind(substring(tmp,first=1,last=first.eq-1),
               substring(tmp,first=first.eq+3))
  tmp <- tmp[tmp[,1]!="",]
  header <- split(tmp[,2],tmp[,1])
  return(header)
}
