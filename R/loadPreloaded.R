#' Load GEO Dataset.
#'
#' \code{loadPreloaded} returns the file with serialized ExpressionSets using
#'     ProtoBuf, that were preloaded on server.
#'
#' @param name String, containing filename. Assuming
#'     that in the directory with preloaded files \code{preloadedDir}
#'     exists file \code{filename.rda} with list of ExpressionSets \code{ess}.
#'
#' @param exactName If you know, that inside file is object with name
#'   \code{exactName}, you can specify it to load only this object.
#'   Otherwise, whole file will be loaded.
#'
#' @return File with ProtoBuf-serialized ExpressionSet-s
#'     that were loaded from specified file.
#'
#' @export
#' @import Biobase
loadPreloaded <- function(name, exactName = NULL) {
  preloadedDir <- getOption("phantasusPreloadedDir")
  if (is.null(preloadedDir)) {
    stop("Specify the directory with presaved files")
  } else if (!dir.exists(preloadedDir)) {
    stop("No such directory")
  }

  fileToLoad <- file.path(preloadedDir, paste0(name, '.rda'))

  if (file.exists(fileToLoad)) {
    x <- load(fileToLoad) # must return the object ess
    loaded <- get(x)

    wrongFormat <- paste("Wrong format.",
                         "File must contain either ExpressionSet or list of ExpressionSets")

    ess <- NULL

    if (class(loaded) == "ExpressionSet") {
      ess <- list()
      ess[[x]] <- loaded
    } else if (class(loaded) != "list") {
      stop(wrongFormat)
    } else {
      ess <- loaded
    }



    files <- list()

    seriesNames <- names(ess)
    if (is.null(seriesNames)) {
      seriesNames <- paste0(name, "_", 1:length(ess))
    }

    if (!is.null(exactName) && !(exactName %in% seriesNames)) {
      stop("There is not such object in this file")
    }

    if (!is.null(exactName)) {
      ess <- list(ess[[exactName]])
    }
    for (i in 1:length(ess)) {
      if (class(ess[[i]]) != "ExpressionSet") {
        stop(wrongFormat)
      }
      assign(paste0("es_", i), ess[[i]], envir = parent.frame())

      files[[seriesNames[[i]]]] <- writeToList(ess[[i]])
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

#' Check names inside preloaded file
#'
#' \code{checkPreloadedNames} checks names of ExpressionSets that
#'   are included in file \code{name}
#'
#' @param name String, containing filename. Assuming
#'     that in the directory with preloaded files \code{preloadedDir}
#'     exists file \code{filename.rda} with list of ExpressionSets \code{ess}.
#'
#' @return Vector of names serialized in JSON format.
#'
#' @export
checkPreloadedNames <- function(name) {
  if (jsonlite::fromJSON(preloadedDirExists())) {
    preloadedDir <- getOption("phantasusPreloadedDir")
    fileToLoad <- file.path(preloadedDir, paste0(name, '.rda'))

    if (file.exists(fileToLoad)) {
      x <- load(fileToLoad) # must return the object ess
      loaded <- get(x)

      answer <- c()

      if (class(loaded) == "ExpressionSet") {
        answer <- c(x)
      } else if (class(loaded) == "list") {
        if (!is.null(names(loaded))) {
          answer <- names(loaded)
        }
        else {
          answer <- paste0(name, "_", 1:length(loaded))
        }
      } else {
        wrongFormat <- paste("Wrong format.",
                             "File must contain either ExpressionSet or list of ExpressionSets")
        stop(wrongFormat)
      }

      jsonlite::toJSON(answer)
    } else {
      stop("No such file")
    }
  } else {
    stop("No such directory")
  }
}

