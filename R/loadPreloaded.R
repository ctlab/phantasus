#' Load GEO Dataset.
#'
#' \code{loadPreloaded} returns the file with serialized ExpressionSets using
#'     ProtoBuf, that were preloaded on server.
#'
#' @param name String, containing filename. Assuming
#'     that in the directory with preloaded files \code{preloadedDir}
#'     exists file \code{filename.rda} with list of ExpressionSets \code{ess}.
#'
#' @return File with ProtoBuf-serialized ExpressionSet-s
#'     that were loaded from specified file.
#'
#' @export
#' @import Biobase
loadPreloaded <- function(name) {
  preloadedDir <- getOption("phantasusPreloadedDir")
  if (is.null(preloadedDir)) {
    stop("Specify the directory with presaved files")
  } else if (!dir.exists(preloadedDir)) {
    stop("No such directory")
  }

  fileToLoad <- file.path(preloadedDir, paste0(name, '.rda'))

  if (file.exists(fileToLoad)) {
    x <- load(fileToLoad) # must return the object ess
    ess <- get(x)

    wrongFormat <- paste("Wrong format.",
                         "File must contain either ExpressionSet or list of ExpressionSets")

    if (class(ess) == "ExpressionSet") {
      ess <- list(ess)
    } else if (class(ess) != "list") {
      stop(wrongFormat)
    }


    files <- list()
    for (i in 1:length(ess)) {
      if (class(ess[[i]]) != "ExpressionSet") {
        stop(wrongFormat)
      }
      assign(paste0("es_", i), ess[[i]], envir = parent.frame())
      seriesName <- paste0(name, '_', i)
      files[[seriesName]] <- writeToList(ess[[i]])
    }
    f <- tempfile(pattern = "gse", tmpdir = getwd(), fileext = ".bin")
    writeBin(protolite::serialize_pb(files), f)
    jsonlite::toJSON(f)
  } else {
    stop("No such file")
  }
}

#' Check existence of phantasusPreloadedDir
#'
#' \code{preloadedDirExists} checks if there is specified
#'   directory with preloaded files.
#'
#' @return Boolean value.
#'
#' @export
preloadedDirExists <- function() {
  preloadedDir <- getOption("phantasusPreloadedDir")
  jsonlite::toJSON(!is.null(preloadedDir) && dir.exists(preloadedDir))
}
