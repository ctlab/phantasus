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
    bin_is_valid <- checkBinValidity(cachedFile, file.info(fileToLoad)$ctime)
    if (!bin_is_valid) {
        writeBin(protolite::serialize_pb(list(layout_version = as.raw(PROTOBUF_LAYOUT_VERSION), ess = files)), cachedFile)
    }

    path <- paste0('/preloaded/', binaryName)

    jsonlite::toJSON(path)

}
#' Generate files for preloaded session from a session link.
#'
#' @param sessionURL String with session link produced by phantasus.
#' @param preloadedName String with name that should be assigned to the session.
#' @param preloadedDir  Path to the directory with preloaded datasets and sessions.
#' @return Function produces two files (\code{preloadedName.rda} with ExpressionSet
#'         and \code{preloadedName.json} with session features) in \code{preloadedDir} folder.
#' @export
#' @examples
#' sessionURL <- "https://ctlab.itmo.ru/phantasus/?session=x063c1b365b9211" # link from 'Get dataset link...' tool in phantasus
#' newName <- "my_session" # user defined name
#' preloadedDir <- "./preloaded" # directory where files will be stored. In order too get access through phantasus web-app should be preloadedDir
#' dir.create(preloadedDir, showWarnings = FALSE)
#' generatePreloadedSession(sessionURL= sessionURL,
#'                          preloadedName = newName,
#'                          preloadedDir = preloadedDir)
#' \dontrun{
#' servePhantasus(preloadedDir=preloadedDir, openInBrowser=FALSE)
#' # open browser manually at http://0.0.0.0:8000/phantasus/index.html?preloaded=my_session
#' }
generatePreloadedSession <- function (sessionURL, preloadedName, preloadedDir) {
    sessionName <-gsub(".+session=([^&]+).*","\\1", sessionURL)
    phantasusPath <- gsub("([^?])(index\\.html)?\\?.*","\\1", sessionURL)
    if (is.null(preloadedDir)) {
        stop("Specify the directory with presaved files")
    } else if (!dir.exists(preloadedDir)) {
        stop("Preloaded directory  does not exist")
    }
    RDataPath = tempfile()
    RDA_req =curl::curl_fetch_disk(url = paste0(phantasusPath,
                                                "ocpu/tmp/", sessionName, "/files/.RData"),
                                   path = RDataPath)

    if (!RDA_req$status_code==200) {
        stop('Invalid session')
    }
    env <- new.env()
    load(RDataPath, envir = env)
    file.remove(RDataPath)

    if (!is.null(env$es)) {
        assign("es", env$es, envir = parent.frame())
        save(es, file = file.path(preloadedDir, paste0(preloadedName, ".rda")))

        if (!is.null(env$heatmapJson)) {
            assign("heatmapJson", env$heatmapJson, envir = parent.frame())
            heatmapPath <- file.path(preloadedDir, paste0(preloadedName,'.json'))
            write(jsonlite::toJSON(heatmapJson, auto_unbox = TRUE, null = 'null'), heatmapPath)
        }
    }
}


