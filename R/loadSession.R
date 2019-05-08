sessionExists <- function(sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- file.path(ocpuRoot, 'ocpu-store', sessionName)

    savedPath <- file.path(sessionPath, 'sess.bin')
    RDataPath <- file.path(sessionPath, '.RData')

    if (file.exists(savedPath)) {
        return (jsonlite::toJSON(list(result=TRUE), auto_unbox = TRUE))
    }

    #-----------------LEGACY-----------------------------------

    if (file.exists(RDataPath)) {
        env <- new.env()
        load(RDataPath, envir = env)
        validForLoading = !is.null(env$es) && !is.null(attr(env$es, 'published'))

        return (jsonlite::toJSON(list(result=validForLoading), auto_unbox = TRUE))
    }

    return (jsonlite::toJSON(list(result=FALSE), auto_unbox = TRUE))
}

publishSession <- function (sessionName, datasetName = sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- file.path(ocpuRoot, 'ocpu-store', sessionName)

    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    if (file.exists(RDataPath)) {
        env <- new.env()
        load(RDataPath, envir = env)

        if (!is.null(env$es)) { # this session is for actual dataset
            result <- list()
            datasetName <- ifelse(nchar(datasetName) > 0, datasetName, sessionName)
            result[[datasetName]] <- writeToList(env$es)
            binaryFile <- file.path(getwd(), 'sess.bin')

            writeBin(protolite::serialize_pb(result), binaryFile)
            assign("es", env$es, envir = parent.frame())
            return (jsonlite::toJSON(basename(binaryFile)))
        }
    }

    stop('Invalid session')
}

loadSession <- function (sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- file.path(ocpuRoot, 'ocpu-store', sessionName)

    savedPath <- paste(sessionPath, 'sess.bin', sep=.Platform$file.sep)
    if (file.exists(savedPath)) {
        return (jsonlite::toJSON(basename(savedPath)))
    }

    #---------------------LEGACY-------------------------------------

    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    if (file.exists(RDataPath)) {
        savedEnv <- load(RDataPath)

        if (!("es" %in% savedEnv) || is.null(attr(es, 'published'))) {
            stop('Invalid session key')
        }

        result <- list()
        result[[sessionName]] <- writeToList(es)
        writeBin(protolite::serialize_pb(result), savedPath)
        return (jsonlite::toJSON(basename(savedPath)))
    }

    stop('Invalid session key')
}
