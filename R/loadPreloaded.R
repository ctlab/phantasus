#' Load GEO Dataset.
#'
#' \code{loadPreloaded} returns the file with serialized ExpressionSets using
#'     ProtoBuf, that were preloaded on server.
#'
#' @param name String, containing filename. Assuming
#'     that in the directory with preloaded files \code{preloadedDir}
#'     exists file \code{filename.rda} with list of ExpressionSets \code{ess}.
#'
#'
#' @return File with ProtoBuf-serialized ExpressionSet-s
#'     that were loaded from specified file.
#'
#' @import Biobase
loadPreloaded <- function(name) {
    preloadedDir <- getOption("phantasusPreloadedDir")
    cacheDir <- getOption("phantasusCacheDir")
    if (is.null(preloadedDir)) {
        stop("Specify the directory with presaved files")
    } else if (!dir.exists(preloadedDir)) {
        stop("No such directory")
    }

    fileToLoad <- file.path(preloadedDir, paste0(name, '.rda'))
    if (!file.exists(fileToLoad)) {
        stop('No such dataset')
    }

    x <- load(fileToLoad) # must return the object ess
    loaded <- get(x)

    wrongFormat <- paste("Wrong format.",
                            "File must contain either ExpressionSet",
                            "or list of ExpressionSets")

    ess <- NULL

    if (class(loaded) == "ExpressionSet") {
        ess <- list()
        ess[[name]] <- loaded
    } else if (class(loaded) == "list") {
        ess <- loaded
        if (!all(unlist(lapply(loaded, class)) == "ExpressionSet")) {
            stop(wrongFormat)
        }
    } else {
        stop(wrongFormat)
    }

    files <- list()
    assign("es", ess[[1]], envir = parent.frame())
    files[[name]] <- writeToList(ess[[1]])
    path <- NULL

    if (is.null(cacheDir)) {
        cacheDir <- tempdir()
    }

    preloadedCacheDir <- file.path(cacheDir, 'preloaded')
    binaryName <- paste0(name, '.bin')
    dir.create(preloadedCacheDir, showWarnings = FALSE, recursive = TRUE)



    cachedFile <- file.path(preloadedCacheDir, binaryName)
    dir.create(dirname(cachedFile), showWarnings = FALSE, recursive = TRUE)
    if (!file.exists(cachedFile)) {
        writeBin(protolite::serialize_pb(files), cachedFile)
    }

    path <- paste0('/preloaded/', binaryName)

    jsonlite::toJSON(path)

}

