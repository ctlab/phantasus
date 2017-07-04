#' @name loadGEO
#' @title Load GEO Dataset
#' @description Function for loading dataset from GEO, creating ExpressionSet from it and exporting it to global environment
#' @param name string, containing GEO identifier of the dataset
#' @param type "GSE" or "GDS"
#' @return file with serialized information on ExpressionSet
#' @export
#' @import Biobase
#' @import GEOquery
loadGEO <- function(name, type = NA) {
  cacheDir <- getOption("phantasusCacheDir")
  if (is.null(cacheDir)) {
    cacheDir <- tempdir()
  } else if (!dir.exists(cacheDir)) {
    dir.create(cacheDir)
  }

  ess <- getES(name, type, destdir = cacheDir)

  writeToList <- function(es) {
    data <- as.matrix(exprs(es)); colnames(data) <- NULL; row.names(data) <- NULL

    pdata <- as.matrix(pData(es)); colnames(pdata) <- NULL; row.names(pdata) <- NULL

    participants <- colnames(es)
    rownames <- rownames(es)

    fdata <- as.matrix(fData(es))
    colnames(fdata) <- NULL
    row.names(fdata) <- NULL

    res <- list(data = data, pdata = pdata,
                fdata = fdata, rownames = rownames,
                colMetaNames = varLabels(phenoData(es)),
                rowMetaNames = varLabels(featureData(es)))
    res
  }


  files <- list()
  for(i in 1:length(ess)) {
    assign(paste("es_", i, sep = ""), ess[[i]], envir = parent.frame())
    seriesName <- paste(name, "-", annotation(ess[[i]]), sep = "")
    files[[seriesName]] <- writeToList(ess[[i]])
  }
  f <- tempfile(pattern = "gse", tmpdir = getwd(), fileext = ".bin")
  writeBin(protolite::serialize_pb(files), f)
  jsonlite::toJSON(f);
}

getGDS <- function(name, destdir = tempdir()) {
  l <- getGEO(name, destdir = destdir)
  table <- slot(l, 'dataTable') # extracting all useful information on dataset
  data <- Table(table)  # extracting table ID_REF | IDENTIFIER/SAMPLE | SAMPLE1 | ...
  columnsMeta <- Columns(table) # phenoData
  sampleNames <- as.vector(columnsMeta[["sample"]])
  rownames <- as.vector(data[["ID_REF"]])
  symbol <- as.vector(data[["IDENTIFIER"]])

  data <- data[,sampleNames] # expression data
  exprs <- as.matrix(data)
  row.names(exprs) <- rownames

  row.names(columnsMeta) <- sampleNames
  # columnsMeta <- columnsMeta[,!(colnames(columnsMeta) %in% c('sample'))]
  pData <- AnnotatedDataFrame(data.frame(columnsMeta, check.names = F))

  fData <- data.frame(matrix(symbol, nrow(exprs), 1));
  colnames(fData) <- "symbol"
  fData <- AnnotatedDataFrame(fData)
  featureNames(fData) <- rownames

  list(ExpressionSet(assayData = exprs, phenoData = pData, featureData = fData))
}

getGSE <- function(name, destdir = tempdir()) {
  GEO <- unlist(strsplit(name, '-'))[1]

  stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  filename <- sprintf("%s_series_matrix.txt.gz", name)
  gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s"

  destfile <- paste(destdir, .Platform$file.sep, filename, sep='')

  infile <- TRUE
  if (!file.exists(destfile)) {
    tryCatch({
        download.file(sprintf(gdsurl, stub, GEO, filename), destfile = destfile)
      }, error = function(e) {
        file.remove(destfile)
      }, warning = function(w) {
        file.remove(destfile)
      }, finally = {
        infile <- file.exists(destfile)
      })
  }

  if (infile) {
    ess <- c(getGEO(filename = destfile, destdir = destdir))
  } else {
    ess <- getGEO(GEO = name, destdir = destdir)
  }

  take <- function(x, n) {
    sapply(x, function(x) { x[[n]] })
  }

  rename <- function(prevName, x) {
    splitted <- strsplit(x, ": ")
    lengths <- sapply(splitted, length)
    if (any(lengths != 2 & lengths != 0)) {
       return(list(name = prevName, x = x))
    }
    splittedFirst <- unique(take(splitted[lengths > 0], 1))
    if (length(splittedFirst) == 1) {

       res = list(name = splittedFirst[1], x = ifelse(lengths == 2,
                                                      take(splitted[lengths == 2], 2),
                                                      NA))

    }
    else {
      res = list(name = prevName, x = x)
    }
    res
  }

  processInputES <- function(es) {
    featureData(es) <- featureData(es)[,grepl("symbol", fvarLabels(es), ignore.case = T)]

    phenoData(es) <- phenoData(es)[,grepl("characteristics", varLabels(es), ignore.case = T)
                                   | (varLabels(es) %in% c("title", "id", "geo_accession"))]

    chr <- varLabels(es)[grepl("characteristics", varLabels(es), ignore.case = T)]

    renamed <- lapply(chr, function(x) { rename(x, as.vector(pData(es)[,x])) })
    phenoData(es) <- phenoData(es)[, !(varLabels(es) %in% chr)]
    pData(es)[,take(renamed,1)] <- take(renamed,2)

    es
  }
  lapply(ess, processInputES)
}


getES <- function(name, type = NA, destdir = tempdir()) {
  if (is.na(type)) {
    type = substr(name, 1, 3)
  }
  possibly.cached <- file.path(destdir, paste(name, '.rda', sep=''))
  if (file.exists(possibly.cached)) {
    load(possibly.cached)
  }
  else {
    if (type == 'GSE') {
      res <- getGSE(name, destdir)
    }
    else if (type == "GDS") {
      res <- getGDS(name, destdir)
    }
    else {
      stop("Incorrect name or type of the dataset")
    }
    if (length(res) > 1) {
      for (i in 1:length(res)) {
        ess <- c(res[[i]])
        save(ess, file = file.path(destdir, paste(name, '-', annotation(res[[i]]), '.rda', sep = '')))
      }
    }
    ess <- res
    save(ess, file = file.path(destdir, paste(name, '.rda', sep = '')))
  }
  ess
}

#' @name checkGPLs
#' @export
checkGPLs <- function(name) {
  stub = gsub("\\d{1,3}$", "nnn", name, perl = TRUE)
  gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
  file.names = GEOquery:::getDirListing(sprintf(gdsurl, stub, name))

  file.names <- file.names[grepl(pattern = paste("^", name, sep = ''), x = file.names)]

  file.names <- unlist(lapply(file.names, function(x) { paste(substr(x, 1, regexpr('_', x) - 1), sep='') }))
  jsonlite::toJSON(file.names)
}

