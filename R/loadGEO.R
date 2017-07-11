#' Load GEO Dataset.
#'
#' \code{loadGEO} returns the file with serialized ExpressionSet using ProtoBuf,
#'   parsed from data downloaded from GEO by identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'   It should start with 'GSE' or 'GDS' and can include exact GPL to annotate dataset,
#'   separated with dash ('-') from the identifier.
#'
#' @param type Type of the dataset: 'GSE' or 'GDS'. If not specified, the function
#'   will take first three letters of \code{name} variable as type.
#'
#' @return File with ProtoBuf-serialized ExpressionSet-s that were downloaded by this identifier.
#'   For GSE-datasets there can be multiple annotations, so in file will be a
#'   list mapping name with GPL to ExpressionSet.
#'
#' @examples
#' loadGEO("GSE27112")
#' loadGEO("GSE27112-GPL1665")
#' loadGEO("GDS4922")
#'
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
                    colMetaNames = varLabels(phenoData(es)),
                    rowMetaNames = varLabels(featureData(es)))
        res
    }


    files <- list()
    for (i in 1:length(ess)) {
        assign(paste0("es_", i), ess[[i]], envir = parent.frame())
        seriesName <- if (!grepl(pattern = "-", name))
            paste0(name, "-", annotation(ess[[i]])) else name
        files[[seriesName]] <- writeToList(ess[[i]])
    }
    f <- tempfile(pattern = "gse", tmpdir = getwd(), fileext = ".bin")
    writeBin(protolite::serialize_pb(files), f)
    jsonlite::toJSON(f)
}

getGDS <- function(name, destdir = tempdir()) {
    l <- getGEO(name, destdir = destdir)

    # extracting all useful information on dataset
    table <- slot(l, "dataTable")

    # extracting table ID_REF | IDENTIFIER/SAMPLE | SAMPLE1 | ...
    data <- Table(table)
    columnsMeta <- Columns(table)  # phenoData
    sampleNames <- as.vector(columnsMeta[["sample"]])
    rownames <- as.vector(data[["ID_REF"]])
    symbol <- as.vector(data[["IDENTIFIER"]])

    data <- data[, sampleNames]  # expression data
    exprs <- as.matrix(data)
    row.names(exprs) <- rownames

    row.names(columnsMeta) <- sampleNames
    pData <- AnnotatedDataFrame(data.frame(columnsMeta, check.names = F))

    fData <- data.frame(matrix(symbol, nrow(exprs), 1))
    colnames(fData) <- "symbol"
    fData <- AnnotatedDataFrame(fData)
    featureNames(fData) <- rownames

    list(ExpressionSet(assayData = exprs,
                       phenoData = pData,
                       featureData = fData))
}

getGSE <- function(name, destdir = tempdir()) {
    GEO <- unlist(strsplit(name, "-"))[1]

    stub <- gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    filename <- sprintf("%s_series_matrix.txt.gz", name)
    gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s"

    destfile <- file.path(destdir, filename)

    infile <- TRUE
    if (!file.exists(destfile)) {
        tryCatch({
            download.file(sprintf(gdsurl, stub, GEO, filename),
                          destfile = destfile)
        }, error = function(e) {
            file.remove(destfile)
        }, warning = function(w) {
            file.remove(destfile)
        }, finally = {
            infile <- file.exists(destfile)
        })
    }

    if (infile) {
        ess <- list(getGEO(filename = destfile, destdir = destdir, AnnotGPL=TRUE))
    } else {
        ess <- getGEO(GEO = name, destdir = destdir, AnnotGPL = TRUE)
    }

    take <- function(x, n) {
        sapply(x, function(x) {
            x[[n]]
        })
    }

    rename <- function(prevName, x) {
        splitted <- strsplit(x, ": ")
        lengths <- sapply(splitted, length)
        if (any(lengths != 2 & lengths != 0)) {
            return(list(name = prevName, x = x))
        }
        splittedFirst <- unique(take(splitted[lengths > 0], 1))
        if (length(splittedFirst) == 1) {

            res <- list(name = splittedFirst[1],
                        x = ifelse(lengths == 2,
                                   take(splitted[lengths == 2], 2),
                                   NA))

        } else {
            res <- list(name = prevName, x = x)
        }
        res
    }

    processInputES <- function(es) {
        featureData(es) <- featureData(es)[, grepl("symbol", fvarLabels(es), ignore.case = T)]

        phenoData(es) <- phenoData(es)[, grepl("characteristics", varLabels(es), ignore.case = T) |
            (varLabels(es) %in% c("title", "id", "geo_accession"))]

        chr <- varLabels(es)[grepl("characteristics", varLabels(es), ignore.case = T)]

        renamed <- lapply(chr, function(x) {
            rename(x, as.vector(pData(es)[, x]))
        })
        phenoData(es) <- phenoData(es)[, !(varLabels(es) %in% chr)]
        pData(es)[, take(renamed, 1)] <- take(renamed, 2)

        es
    }
    lapply(ess, processInputES)
}


getES <- function(name, type = NA, destdir = tempdir()) {
    if (is.na(type)) {
        type = substr(name, 1, 3)
    }
    possibly.cached <- file.path(destdir, paste0(name, ".rda"))
    if (file.exists(possibly.cached)) {
        load(possibly.cached)
    } else {
        if (type == "GSE") {
            res <- getGSE(name, destdir)
        } else if (type == "GDS") {
            res <- getGDS(name, destdir)
        } else {
            stop("Incorrect name or type of the dataset")
        }
        if (length(res) > 1) {
            for (i in 1:length(res)) {
                ess <- c(res[[i]])
                save(ess, file = file.path(destdir, paste0(name, "-", annotation(res[[i]]), ".rda")))
            }
        }
        ess <- res
        save(ess, file = file.path(destdir, paste0(name, ".rda")))
    }
    ess
}

#' Check possible annotations for GEO Dataset.
#'
#' \code{checkGPLs} returns GPL-names for
#'   the specified GEO identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'
#' @return Vector of filenames serialized in JSON format.
#'   If there is only one GPL for that dataset, the function will
#'   return \code{name}.
#'
#' @examples
#' checkGPLs('GSE27112')
#' checkGPLs('GSE14308')
#'
#' @export
checkGPLs <- function(name) {
    if (substr(name, 1, 3) == "GDS") {
      return(jsonlite::toJSON(name))
    }
    stub <- gsub("\\d{1,3}$", "nnn", name, perl = TRUE)
    gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
    file.names <- GEOquery:::getDirListing(sprintf(gdsurl, stub, name))

    file.names <- file.names[grepl(pattern = paste0("^", name), x = file.names)]

    file.names <- unlist(lapply(file.names, function(x) {
        paste0(substr(x, 1, regexpr("_", x) - 1))
    }))
    jsonlite::toJSON(file.names)
}

