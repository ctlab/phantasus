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
#'
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
    for (i in seq_along(ess)) {
        assign(paste0("es_", i), ess[[i]], envir = parent.frame())
        seriesName <- if (!grepl(pattern = "-", name) && length(ess) > 1)
            paste0(name, "-", annotation(ess[[i]])) else name
        files[[seriesName]] <- writeToList(ess[[i]])
    }
    f <- tempfile(pattern = "gse", tmpdir = getwd(), fileext = ".bin")
    writeBin(protolite::serialize_pb(files), f)
    jsonlite::toJSON(f)
}

#' Load ExpressionSet from GEO Datasets
#'
#'\code{getGDS} return the ExpressionSet object corresponding
#'     to GEO Dataset identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'     It should start with 'GSE' or 'GDS' and can include exact GPL
#'     to annotate dataset, separated with dash ('-') from the identifier.
#'
#' @param destdir Directory for caching loaded Series and GPL
#'     files from GEO database.
#'
#' @param mirrorPath URL string which specifies the source of matrices.
#'
#' @return ExpressionSet object wrapped in list, that was available by given
#'     in \code{name} variable GEO identifier.
#'
#' @examples
#' getGDS('GDS4922')
#'
#' @export
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

loadFromARCHS4 <- function(es, archs4_gpls, destdir) {
    if (!es@annotation %in% archs4_gpls$gpl) {
        return(es)
    }

    destfile <- file.path(destdir, archs4_gpls$file[match(es@annotation, archs4_gpls$gpl)])
    if (!file.exists(destfile)) {
        warning(paste0(
            "Found GPL supported by ARCHS4 but not corresponding file available at ",
            destfile))
        return(es)
    }

    samples <- h5read(destfile, "meta/Sample_geo_accession")
    genes <- as.character(h5read(destfile, "meta/genes"))

    sampleIndexes <- match(es$geo_accession,
                           samples)

    expression <- h5read(destfile,
                         "data/expression",
                         index=list(seq_along(genes),
                                    stats::na.omit(sampleIndexes)))
    rownames(expression) <- genes
    colnames(expression) <- colnames(es)[!is.na(sampleIndexes)]
    H5close()

    es2 <- ExpressionSet(assayData = expression,
                         phenoData = phenoData(es[, !is.na(sampleIndexes)]),
                         annotation = annotation(es))
    fData(es2) <- cbind(fData(es2), "Gene symbol"=rownames(es2))
    return(es2)
}

filterFeatureAnnotations <- function(es) {
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
    } else if ("ID" %in% fvarLabels(es)) {
        fvarsToKeep <- c(fvarsToKeep, "ID")
    } else {
        fvarsToKeep <- c(fvarsToKeep, grep("entrez",
                                           fvarLabels(es),
                                           ignore.case = TRUE,
                                           value = TRUE))
    }

    featureData(es) <- featureData(es)[, fvarsToKeep]

    es
}

filterPhenoAnnotations <- function(es) {
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

    parsePData <- function(old.phenodata) {
        old.pdata <- pData(old.phenodata)
        labels <- varLabels(old.phenodata)

        new.pdata <- as.data.frame(matrix(NA, nrow = nrow(old.pdata), ncol = 0))


        for (i in seq_len(ncol(old.pdata))) {
            splitted <- strsplit(as.vector(old.pdata[[i]]), ':')
            lengths <- sapply(splitted, length)
            if (any(lengths != 2 & lengths != 0)) {
                new.pdata[[labels[i]]] <- old.pdata[[i]]
            } else {
                zeros <- which(lengths == 0)

                splitted[zeros] <- replicate(length(zeros), list(c(NA, NA)))

                newnames <- unique(trimws(take(splitted, 1)))
                newnames <- newnames[which(!is.na(newnames))]

                for (j in seq_along(newnames)) {
                    name <- newnames[j]
                    if (!(name %in% names(new.pdata))) {
                        new.pdata[[name]] <- replicate(nrow(new.pdata), NA)
                    }
                    indices <- which(name == trimws(take(splitted, 1)))
                    new.pdata[[name]][indices] <- trimws(take(splitted, 2)[indices])
                }
            }
        }
        AnnotatedDataFrame(new.pdata)
    }

    if (ncol(es) > 0) {
        phenoData(es) <- parsePData(phenoData(es))
    }

    es
}


#' Load ExpressionSet from GEO Series
#'
#'\code{getGSE} return the ExpressionSet object(s) corresponding
#'     to GEO Series Identifier.
#'
#' @param name String, containing GEO identifier of the dataset.
#'     It should start with 'GSE' or 'GDS' and can include exact GPL
#'     to annotate dataset, separated with dash ('-') from the identifier.
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
#'     getGSE('GSE14308', destdir = 'cache')
#'     getGSE('GSE27112')
#' }
#' getGSE('GSE53986')
#'
#' @export
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
        if (length(gpls) == 0) {
            stop(paste("Dataset", name, "not found"))
        }
        ess <- list()
        for (i in 1:length(gpls)) {
          ess[[gpls[[i]]]] <- getGSE(gpls[[i]], destdir = destdir, mirrorPath = mirrorPath)[[1]]
        }
        return(ess)
    }


    archs4_gpls <- rbind(
        data.frame(gpl=c("GPL13112", "GPL17021", "GPL15103"),
                   file="mouse_matrix.h5"),
        data.frame(gpl=c("GPL11154", "GPL16791", "GPL10999", "GPL9115", "GPL18460"),
                   file="human_matrix.h5"))

    ess <- lapply(ess, loadFromARCHS4, archs4_gpls=archs4_gpls, destdir=destdir)

    ess <- lapply(ess, filterFeatureAnnotations)

    ess <- lapply(ess, filterPhenoAnnotations)

    ess
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
#'     getES('GSE14308', type = 'GSE', destdir = 'cache')
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


listCachedESs <- function(destdir) {
    res <- list.files(destdir, pattern=".*\\.rda")
    res <- grep("\\.gz\\.rda$", res, invert = TRUE, value = TRUE)
    res <- grep("^(GSE|GDS)", res, value = TRUE)
    res <- sub("\\.rda$", "", res)
    res
}

#' Reparse cached expression sets from GEO.
#'
#' The function should be used on phantasus version updates that change
#' behavior of loading datasets from GEO. It finds all the datasets
#' that were cached and runs `getES` for them again. The function
#' uses cached Series and other files from GEO.
#'
#' @param destdir Directory used for caching loaded Series files from GEO database.
#'
#' @param mirrorPath URL string which specifies the source of matrices.
#'
#' @return vector of previously cached GSE IDs
#'
#' @examples
#' reparseCachedESs(destdir=tempdir())
#'
#' @export
reparseCachedESs <- function(destdir,
                                mirrorPath = "https://ftp.ncbi.nlm.nih.gov") {
    toReparse <- listCachedESs(destdir)

    for (name in toReparse) {
        message(paste0("Reparsing dataset ", name))
        destfile <- file.path(destdir, paste0(name, ".rda"))
        bakfile <- paste0(destfile, ".bak")
        tryCatch({
            file.rename(destfile, bakfile)
            getES(name, destdir = destdir, mirrorPath = mirrorPath)
            file.remove(bakfile)
        }, error = function(e) {
            message(paste0("Error occured while reparsing, old file stored as ",
                           bakfile))
        })
    }
    return(toReparse)
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
#' \dontrun{
#' checkGPLs('GSE27112')
#' checkGPLs('GSE14308')
#' }
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
                if (length(file.names) == 1) {
                    file.names <- c(name)
                }
                gpls <- file.names
            }

            return(jsonlite::toJSON(gpls))
        }
    },
    error = function(e) {
        message(paste("Problems establishing connection.",
                        "Trying to find corresponding files in cache."))

        files <- list.files(path = cacheDir)

        corresponding <- take(sapply(files[grep(x = files,
                                    pattern = paste0(name, "[-_].*(gz|rda)$"))],
                                    FUN = function(x) { strsplit(x, ".", fixed = TRUE) }), 1)
        gpls <- unique(take(sapply(corresponding,
                            FUN = function(x) { strsplit(x, "_") }), 1))
        if (length(gpls) == 0) {
            warning("No corresponding files were found")
        }
        return(jsonlite::toJSON(gpls))
    })
}
