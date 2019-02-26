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

    geoDir <- getGEODir(name, cacheDir)
    binaryName <- paste0(name, '.bin')
    filePath <- file.path(geoDir, binaryName)
    urls <- c()


    ess <- getES(name, type, destdir = cacheDir, mirrorPath = mirrorPath)

    files <- list()
    assign("es", tail(ess, n=1)[[1]], envir = parent.frame())
    if (!file.exists(filePath)) {
        for (i in seq_along(ess)) {
            seriesName <- if (!grepl(pattern = "-", name) && length(ess) > 1)
                paste0(name, "-", annotation(ess[[i]])) else name
            files[[seriesName]] <- writeToList(ess[[i]])
        }
        writeBin(protolite::serialize_pb(files), filePath)
    }

    urls <-c(urls, unlist(strsplit(filePath, cacheDir))[2])
    jsonlite::toJSON(urls)
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

    fullGEODirPath <- getGEODir(name, destdir)
    destfile <- file.path(fullGEODirPath, filename)

    infile <- file.exists(destfile)

    if (!infile) {
        tempDestFile <- tempfile(paste0(filename, ".load"), tmpdir=fullGEODirPath)
        tryCatch({
            utils::download.file(sprintf(gdsurl, mirrorPath,
                                         stub, name, filename),
                            destfile = tempDestFile,
                            method="libcurl")
            file.rename(tempDestFile, destfile)
            infile <- TRUE
        },
        error = function(e) {
            file.remove(tempDestFile)
        },
        warning = function(w) {
            file.remove(tempDestFile)
        })
    } else {
        message(paste("Loading from locally found file", destfile))
    }

    l <- suppressWarnings(getGEO(GEO = name, destdir = fullGEODirPath, AnnotGPL = TRUE))

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

    fData <- data.frame(id=rownames, symbol=symbol, row.names = rownames, stringsAsFactors = FALSE)
    fData <- AnnotatedDataFrame(fData)

    mabstract=ifelse(is.null(Meta(l)$description),"",Meta(l)$description)
    mpubmedids=ifelse(is.null(Meta(l)$pubmed_id),"",Meta(l)$pubmed_id)
    mtitle=ifelse(is.null(Meta(l)$title),"",Meta(l)$title)

    list(ExpressionSet(assayData = exprs,
                       phenoData = pData,
                       featureData = fData,
                       experimentData=new("MIAME",
                                          abstract=mabstract,
                                          title=mtitle,
                                          pubMedIds=mpubmedids,
                                          other=Meta(l))))
}


#' Returns list of ARCHS4 hdf5 files with expression data
#' @param cacheDir base directory for cache
#' @return list of .h5 files
getArchs4Files <- function(cacheDir) {
    list.files(paste(file.path(cacheDir), 'archs4', sep = .Platform$file.sep), '*.h5', full.names = TRUE)
}


#' Loads expression data from ARCHS4 count files.
#' Only sapmles with counted expression are kept.
#' If es already containts expression data it is returned as is.
#' @param es ExpressionSet from GEO to check for expression in ARCHS4
#' @param archs4_files list of available .h5 files from ARCHS4 project
#' @return either original es or an ExpressionSet with loaded count data from ARCHS4
loadFromARCHS4 <- function(es, archs4_files) {
    if (nrow(es) > 0 ) {
        return(es)
    }
    for (destfile in archs4_files) {

        samples <- h5read(destfile, "meta/Sample_geo_accession")
        sampleIndexes <- match(es$geo_accession,
                               samples)

        if (sum(!is.na(sampleIndexes)) == 0) {
            # no needed samples in this file
            H5close()
            next
        }

        genes <- as.character(h5read(destfile, "meta/genes"))

        expression <- h5read(destfile,
                             "data/expression",
                             index=list(seq_len(length(genes)),
                                        stats::na.omit(sampleIndexes)))
        rownames(expression) <- genes
        colnames(expression) <- colnames(es)[!is.na(sampleIndexes)]
        H5close()

        es2 <- ExpressionSet(assayData = expression,
                             phenoData = phenoData(es[, !is.na(sampleIndexes)]),
                             annotation = annotation(es),
                             experimentData = experimentData(es)
                            )
        fData(es2) <- cbind(fData(es2), "Gene symbol"=rownames(es2))
        return(es2)
    }

    return(es)
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

    if (!any(sapply(fData(es),
                    function(x) identical(rownames(es), as.character(x))
                    ))) {
        fData(es) <- cbind("id"=rownames(es), fData(es))
    }

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
        rownames(new.pdata) <- rownames(old.pdata)
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
    gseurl <- "%s/geo/series/%s/%s/matrix/%s"

    fullGEODirPath <- getGEODir(name, destdir)
    destfile <- file.path(fullGEODirPath, filename)

    infile <- file.exists(destfile)

    if (!infile) {
        tempDestFile <- tempfile(paste0(filename, ".load"), tmpdir=fullGEODirPath)
        tryCatch({
            utils::download.file(sprintf(gseurl, mirrorPath,
                                         stub, GEO, filename),
                                    destfile = tempDestFile,
                                 method="libcurl")
            file.rename(tempDestFile, destfile)
            infile <- TRUE
        },
        error = function(e) {
            file.remove(tempDestFile)
        },
        warning = function(w) {
            file.remove(tempDestFile)
        })
    } else {
        message(paste("Loading from locally found file", destfile))
    }

    if (infile && file.size(destfile) > 0) {
        ess <- list(suppressWarnings(getGEO(filename = destfile,destdir = fullGEODirPath, getGPL = FALSE, AnnotGPL = FALSE)))
        for (i in seq_len(length(ess))) {
            ess[[i]] <- annotateFeatureData(ess[[i]], destdir)
        }
    } else {
        gpls <- fromJSON(checkGPLs(name))
        if (length(gpls) == 0) {
            stop(paste("Dataset", name, "not found"))
        }
        if (length(gpls) == 1 && gpls == name) {
            stop(paste("Can't download dataset ", name))

        }
        ess <- list()
        for (i in 1:length(gpls)) {
          ess[[gpls[[i]]]] <- getGSE(gpls[[i]], destdir = destdir, mirrorPath = mirrorPath)[[1]]
        }
        return(ess)
    }


    archs4_files <- getArchs4Files(destdir)
    if (length(archs4_files) > 0)  {
        ess <- lapply(ess, loadFromARCHS4, archs4_files=archs4_files)
    }

    ess <- lapply(ess, filterFeatureAnnotations)

    ess <- lapply(ess, filterPhenoAnnotations)

    ess <- lapply(ess, inferCondition)

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
    geoDir <- getGEODir(name, destdir)
    possibly.cached <- file.path(geoDir, paste0(name, ".rda"))
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
                destfile <- file.path(geoDir,
                                      paste0(name,
                                             "-",
                                             annotation(res[[i]]),
                                             ".rda"))
                message(paste("Cached dataset to ", destfile))
                save(ess, file = destfile)
            }
        }
        ess <- res
        destfile <- file.path(geoDir, paste0(name, ".rda"))
        message(paste("Cached dataset to ", destfile))
        save(ess, file = destfile)
    }
    ess
}


listCachedESs <- function(destdir) {
    correctFileName <-

    res <- list.files(destdir, pattern=".*\\.rda$", recursive = TRUE)
    res <- grep("\\.gz\\.rda$", res, invert = TRUE, value = TRUE)
    res <- res[grep("^(GSE|GDS)", basename(res), value = FALSE)]
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

    for (path in toReparse) {
        name <- basename(path)
        geoDir <- getGEODir(name, destdir)
        message(paste0("Reparsing dataset ", name))
        destfile <- file.path(geoDir, paste0(name, ".rda"))
        bakfile <- paste0(destfile, ".bak")
        tryCatch({
            file.rename(destfile, bakfile)
            getES(name, destdir = destdir, mirrorPath = mirrorPath)
            file.remove(bakfile)
            file.remove(file.path(geoDir, paste0(name, ".bin")))
        }, error = function(e) {
            message(paste0("Error occured while reparsing, old file stored as ",
                           bakfile))
        })
    }
    return(basename(toReparse))
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

removeRepeatWords <- function(titles) {
    titles_without_repeat_words <- titles
    repeat_words <- regmatches(titles, regexpr("(?![-+}\\]\\)])(\\W|_)*\\w*$", titles, ignore.case = TRUE, perl = TRUE))
    lsuff <- lcSuffix(repeat_words, ignore.case = TRUE)
    if ((all(stringr::str_length(lsuff) == stringr::str_length(repeat_words))) & (all(sub("(\\W*|_)*\\w*$", "", titles, ignore.case = TRUE, perl = TRUE) != ""))) {
        titles_without_repeat_words <- sub("(?![-+}\\]\\)])(\\W|_)*\\w*$", "", titles, ignore.case = TRUE, perl = TRUE)
    }
    if (all(grepl("-$", titles_without_repeat_words, ignore.case = TRUE))) titles_without_repeat_words <- sub("-$", "", titles_without_repeat_words, ignore.case = TRUE, perl = TRUE)
    return(titles_without_repeat_words)
}

inferConditionImpl <- function(gse_titles) {
    inferCondition <- gse_titles
    rep_num <- NULL
    if ((length(inferCondition) > 40) | (length(inferCondition) < 3))
    {
        return(list())
    } else if (! all(grepl("[A-z]", inferCondition, ignore.case = TRUE)))
    {
        return(list())
    } else if (all(grepl(" vs[. ]{1}", inferCondition)))
    {
        return(list())
    } else
    {
        lsuff <- lcSuffix(inferCondition)
        suff_length <- stringr::str_length(lsuff)
        if (suff_length > 1) inferCondition <- stringr::str_sub(inferCondition, 1, -suff_length-1)
        if (all(grepl("((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+", inferCondition, ignore.case = TRUE)))
        {
            sample <- regmatches(inferCondition, regexpr("((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+", inferCondition, ignore.case = TRUE))
            rep_num <- regmatches(sample, regexpr("\\d+$", sample))
            inferCondition <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*\\d+(\\W|_)*", "", inferCondition, ignore.case = TRUE, perl = TRUE)

        }
        else if (all(grepl("((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*[A-z]$", inferCondition, ignore.case = TRUE)))
        {
            sample <- regmatches(inferCondition, regexpr("((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*[A-z]$", inferCondition, ignore.case = TRUE))
            rep_num <- regmatches(sample, regexpr("[A-z]$", sample))
            inferCondition <- sub("(?![-+}\\]\\)])(\\W|_)*((bio[A-z]*)|(tech[A-z]*))?[ .#_-]*((rep.*)|(set.*)|(sample)|(exp.*)|(case))[ .#_-]*[A-z]$", "", inferCondition, ignore.case = TRUE, perl = TRUE)

        }
        else if (length(unique(sub("[- \\.#_]*\\d+$", "", inferCondition))) == 1) {
            return(list())
        }
        else if (all(grepl("[- \\.#_]*\\d+$", inferCondition)) & (length(unique(regmatches(inferCondition, regexpr("\\d+$", inferCondition)))) > 1))
        {
            sample <- sub("\\d+$", "", inferCondition)
            if (length(unique(regmatches(sample, regexpr("[- \\.#_]$", sample)))) > 1)
            {
                return(list())
            }
            else
            {
                rep_num <- regmatches(inferCondition, regexpr("\\d+$", inferCondition))
                inferCondition <- removeRepeatWords(sub("[ \\.#_]*\\d+$", "", inferCondition))

            }
        }
        else return(list())
    }

    if ((length(rep_num) > 1) & (length(unique(inferCondition)) < length(unique(gse_titles))) & (length(unique(inferCondition)) > 1))
        return(list(condition=inferCondition, replicate=rep_num))
    else return(list())
}

inferCondition <- function(es) {
    newAnnot <- inferConditionImpl(es$title)
    if (length(newAnnot) == 2) {
        pData(es)$condition <- newAnnot$condition
        pData(es)$replicate <- newAnnot$replicate
    }
    es
}


getGPLAnnotation <- function (GPL, destdir = tempdir()) {
    GPL <- toupper(GPL)
    stub = gsub('\\d{1,3}$','nnn',GPL,perl=TRUE)
    GPLDirPath <- '%s/geo/platforms/%s/%s/annot'
    fullGPLDirPath <- file.path(sprintf(GPLDirPath, destdir, stub, GPL))

    dir.create(fullGPLDirPath, showWarnings = FALSE, recursive = TRUE)

    return(suppressWarnings(getGEO(GPL, destdir = fullGPLDirPath, AnnotGPL = TRUE)))
}

annotateFeatureData <- function (es, destdir = tempdir()) {
    platform <- levels(es$platform_id)[1]
    platformParsed <- getGPLAnnotation(platform, destdir)

    #https://github.com/seandavi/GEOquery/blob/master/R/parseGEO.R#L569
    #############################
    vmd <- Columns(platformParsed)
    dat <- Table(platformParsed)
    ## GEO uses "TAG" instead of "ID" for SAGE GSE/GPL entries, but it is apparently
    ##     always the first column, so use dat[,1] instead of dat$ID
    ## The next line deals with the empty GSE
    tmpnames=character(0)
    if(ncol(dat)>0) {
        tmpnames=as.character(dat[,1])
    }
    ## Fixed bug caused by an ID being "NA" in GSE15197, for example
    tmpnames[is.na(tmpnames)]="NA"
    rownames(dat) <- make.unique(tmpnames)
    ## Apparently, NCBI GEO uses case-insensitive matching
    ## between platform IDs and series ID Refs ???
    dat <- dat[match(tolower(rownames(es)),tolower(rownames(dat))),]
    # Fix possibility of duplicate column names in the
    # GPL files; this is prevalent in the Annotation GPLs
    rownames(vmd) <- make.unique(colnames(Table(platformParsed)))
    colnames(dat) <- rownames(vmd)
    ##############################

    featureData(es) <- new('AnnotatedDataFrame',data=dat,varMetadata=vmd)
    es
}

getGEODir <- function (name, destdir) {
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
