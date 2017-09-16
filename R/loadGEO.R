#' Load GEO Dataset.
#'
#' \code{loadGEO} returns the file with serialized ExpressionSet using
#'     ProtoBuf, parsed from data downloaded from GEO by identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'     It should start with 'GSE' or 'GDS' and can include exact GPL
#'     to annotate dataset, separated with dash ('-') from the identifier.
#'
#' @param type Type of the dataset: 'GSE' or 'GDS'. If not specified,
#'     the function will take first three letters
#'     of \code{name} variable as type.
#'
#' @return File with ProtoBuf-serialized ExpressionSet-s
#'     that were downloaded by this identifier.
#'     For GSE-datasets there can be multiple annotations, so in file will be a
#'     list mapping name with GPL to ExpressionSet.
#'
#' @examples
#' \dontrun{
#'     loadGEO("GSE27112")
#'     loadGEO("GDS4922")
#' }
#' loadGEO("GSE27112-GPL6885")
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

    mirrorPath <- getOption('phantasusMirrorPath')
    if (is.null(mirrorPath)) {
      mirrorPath <- "https://ftp.ncbi.nlm.nih.gov"
    }

    ess <- getES(name, type, destdir = cacheDir, mirrorPath = mirrorPath)

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

getGDS <- function(name, destdir = tempdir(),
                   mirrorPath = "https://ftp.ncbi.nlm.nih.gov") {
    stub <- gsub("\\d{1,3}$", "nnn", name, perl = TRUE)
    filename <- sprintf("%s.soft.gz", name)
    gdsurl <- "%s/geo/datasets/%s/%s/soft/%s"

    destfile <- file.path(destdir, filename)

    infile <- FALSE
    if (!file.exists(destfile)) {
        tryCatch({
            utils::download.file(sprintf(gdsurl, mirrorPath,
                                         stub, name, filename),
                            destfile = destfile)
            infile <- TRUE
        },
        error = function(e) {
            file.remove(destfile)
        },
        warning = function(w) {
            file.remove(destfile)
        })
    } else {
        message(paste("Loading from locally found file", destfile))
    }

    if (infile) {
        l <- suppressWarnings(getGEO(filename = destfile,
                                        destdir = destdir,
                                        AnnotGPL = TRUE))
    } else {
        l <- suppressWarnings(getGEO(GEO = name,
                                    destdir = destdir,
                                    AnnotGPL = TRUE))
    }

    # extracting all useful information on dataset
    table <- methods::slot(l, "dataTable")

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
    pData <- AnnotatedDataFrame(data.frame(columnsMeta, check.names = FALSE))

    fData <- data.frame(matrix(symbol, nrow(exprs), 1))
    colnames(fData) <- "symbol"
    fData <- AnnotatedDataFrame(fData)
    featureNames(fData) <- rownames

    list(ExpressionSet(assayData = exprs,
                        phenoData = pData,
                        featureData = fData))
}

#' @import rhdf5
getGSE <- function(name, destdir = tempdir(),
                   mirrorPath = "https://ftp.ncbi.nlm.nih.gov") {
    GEO <- unlist(strsplit(name, "-"))[1]

    stub <- gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    filename <- sprintf("%s_series_matrix.txt.gz", name)
    gdsurl <- "%s/geo/series/%s/%s/matrix/%s"

    destfile <- file.path(destdir, filename)

    infile <- file.exists(destfile)
    if (!file.exists(destfile)) {
        tryCatch({
            utils::download.file(sprintf(gdsurl, mirrorPath,
                                         stub, GEO, filename),
                                    destfile = destfile)
            infile <- TRUE
        },
        error = function(e) {
            file.remove(destfile)
        },
        warning = function(w) {
            file.remove(destfile)
        })
    } else {
        message(paste("Loading from locally found file", destfile))
    }

    if (infile) {
        ess <- list(suppressWarnings(getGEO(filename = destfile,
                                                destdir = destdir,
                                                AnnotGPL = TRUE)))
    } else {
        gpls <- fromJSON(checkGPLs(name))
        ess <- list()
        for (i in 1:length(gpls)) {
          ess[[gpls[[i]]]] <- getGSE(gpls[[i]], destdir = destdir, mirrorPath = mirrorPath)[[1]]
        }
    }


    archs4_gpls <- rbind(
        data.frame(gpl=c("GPL13112", "GPL17021", "GPL15103"),
                   file="mouse_matrix.h5"),
        data.frame(gpl=c("GPL11154", "GPL16791", "GPL10999", "GPL9115", "GPL18460"),
                   file="human_matrix.h5"))

    for (i in seq_along(ess)) {
        es <- ess[[i]]
        if (es@annotation %in% archs4_gpls$gpl) {
            destfile <- file.path(destdir, archs4_gpls$file[match(es@annotation, archs4_gpls$gpl)])
            if (!file.exists(destfile)) {
                warning(paste0(
                    "Found GPL supported by ARCHS4 but not corresponding file available at ",
                    destfile))
                next
            }

            samples <- h5read(destfile, "meta/Sample_geo_accession")
            genes <- as.character(h5read(destfile, "meta/genes"))

            sampleIndexes <- match(es$geo_accession,
                                   samples)

            expression <- h5read(destfile,
                                 "data/expression",
                                 index=list(seq_along(genes),
                                            na.omit(sampleIndexes)))
            rownames(expression) <- genes
            colnames(expression) <- colnames(es)[!is.na(sampleIndexes)]
            H5close()

            es2 <- ExpressionSet(assayData = expression,
                                 phenoData = phenoData(es[, !is.na(sampleIndexes)]),
                                 annotation = annotation(es))
            fData(es2) <- cbind(fData(es2), "Gene symbol"=rownames(es2))
            ess[[i]] <- es2

        }
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
        fvarsToKeep <- c()
        if ("Gene symbol" %in% fvarLabels(es)) {
            fvarsToKeep <- c(fvarsToKeep, "Gene symbol")
        } else {
            fvarsToKeep <- c(fvarsToKeep, grep("symbol",
                                                fvarLabels(es),
                                                ignore.case = TRUE,
                                                value = TRUE))
        }

        if ("Gene ID" %in% fvarLabels(es)) {
            fvarsToKeep <- c(fvarsToKeep, "Gene ID")
        } else  {
            fvarsToKeep <- c(fvarsToKeep, grep("entrez",
                                                fvarLabels(es),
                                                ignore.case = TRUE,
                                                value = TRUE))
        }

        featureData(es) <- featureData(es)[, fvarsToKeep]

        phenoData(es) <- phenoData(es)[,
                                        grepl("characteristics",
                                                varLabels(es),
                                                ignore.case = TRUE) |
                                        (varLabels(es) %in% c("title",
                                                                "id",
                                                                "geo_accession"
                                                                ))]

        chr <- varLabels(es)[grepl("characteristics",
                                    varLabels(es),
                                    ignore.case = TRUE)]

        renamed <- lapply(chr, function(x) {
            rename(x, as.vector(pData(es)[, x]))
        })
        phenoData(es) <- phenoData(es)[, !(varLabels(es) %in% chr)]
        pData(es)[, take(renamed, 1)] <- take(renamed, 2)

        es
    }
    lapply(ess, processInputES)
}

#' Load ExpressionSet by GEO identifier
#'
#'\code{getES} return the ExpressionSet object(s) corresponding
#'     to GEO identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'     It should start with 'GSE' or 'GDS' and can include exact GPL
#'     to annotate dataset, separated with dash ('-') from the identifier.
#'
#' @param type Type of the dataset: 'GSE' or 'GDS'. If not specified,
#'     the function will take first three letters
#'     of \code{name} variable as type.
#'
#' @param destdir Directory for caching loaded Series and GPL
#'     files from GEO database.
#'
#' @param mirrorPath URL string which specifies the source of matrices.
#'
#' @return List of ExpressionSet objects, that were available by given
#'     in \code{name} variable GEO identifier.
#'
#' @examples
#' \dontrun{
#'     getES('GSE14308', type = 'GSE', destdir = file.path(getwd(), 'cache'))
#'     getES('GSE27112')
#' }
#' getES('GDS4922')
#'
#' @export
getES <- function(name, type = NA, destdir = tempdir(),
                  mirrorPath = "https://ftp.ncbi.nlm.nih.gov") {
    if (is.na(type)) {
        type <- substr(name, 1, 3)
    }
    possibly.cached <- file.path(destdir, paste0(name, ".rda"))
    if (file.exists(possibly.cached)) {
        load(possibly.cached)
        message(paste("Loaded from locally cached parsed file", possibly.cached))
    } else {
        if (type == "GSE") {
            res <- getGSE(name, destdir, mirrorPath)
        } else if (type == "GDS") {
            res <- getGDS(name, destdir, mirrorPath)
        } else {
            stop("Incorrect name or type of the dataset")
        }
        if (length(res) > 1) {
            for (i in 1:length(res)) {
                ess <- list(res[[i]])
                save(ess, file = file.path(destdir,
                                            paste0(name,
                                                    "-",
                                                    annotation(res[[i]]),
                                                    ".rda")))
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
#'     the specified GEO identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'
#' @return Vector of filenames serialized in JSON format.
#'     If there is only one GPL for that dataset, the function will
#'     return \code{name}.
#'
#' @examples
#' checkGPLs('GSE27112')
#' checkGPLs('GSE14308')
#'
#' @export
checkGPLs <- function(name) {
    mirrorPath <- getOption('phantasusMirrorPath')
    if (is.null(mirrorPath)) {
      mirrorPath <- "https://ftp.ncbi.nlm.nih.gov"
    }

    cacheDir <- getOption("phantasusCacheDir")

    if (is.null(cacheDir)) {
      cacheDir <- tempdir()
    } else if (!dir.exists(cacheDir)) {
      dir.create(cacheDir)
    }

    type <- substr(name, 1, 3)
    assertthat::assert_that( (type == "GDS" || type == "GSE")
                            && nchar(name) >= 4)

    stub <- gsub("\\d{1,3}$", "nnn", name, perl = TRUE)
    gdsurl <- "%s/geo/%s/%s/%s/"

    url <- sprintf(gdsurl, mirrorPath,
                   if (type == "GDS") "datasets" else "series", stub, name)

    gpls <- c()

    tryCatch({
        httr::GET(url)
        if (httr::status_code(httr::GET(url)) == 404) {
            warning("No such dataset")
            return(jsonlite::toJSON(c()))
        } else {
            if (type == "GDS") {
                gpls <- c(name)
            } else {
                file.names <- GEOquery:::getDirListing(paste0(url, "matrix/"))

                file.names <- file.names[grepl(pattern = paste0("^", name),
                                            x = file.names)]

                file.names <- unlist(lapply(file.names, function(x) {
                    paste0(substr(x, 1, regexpr("_", x) - 1))
                }))
                gpls <- file.names
            }
        }
    },
    error = function(e) {
        message(paste("Problems establishing connection.",
                        "Trying to find corresponding files in cache."))

        files <- list.files(path = cacheDir)

        corresponding <- files[grep(x = files,
                                    pattern = paste0(name, "[-_].*gz$"))]
        gpls <- take(sapply(corresponding,
                            FUN = function(x) { strsplit(x, "_") }), 1)
        if (length(gpls) == 0) {
            warning("No corresponding files were found")
        }
    })

    return(jsonlite::toJSON(gpls))

}
