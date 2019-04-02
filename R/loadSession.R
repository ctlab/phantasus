sessionExists <- function(sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)

    savedPath <- paste(sessionPath, 'sess.bin', sep=.Platform$file.sep)
    if (file.exists(savedPath)) {
        return (jsonlite::toJSON(list(result=TRUE), auto_unbox = TRUE))
    }

    return (jsonlite::toJSON(list(result=FALSE), auto_unbox = TRUE))
}

publishSession <- function (sessionName, datasetName = sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)

    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    if (file.exists(RDataPath)) {
        env <- new.env()
        load(RDataPath, envir = env)

        if (!is.null(env$es)) { # this session is for actual dataset
            result <- list()
            datasetName <- ifelse(!is.null(datasetName) && nchar(datasetName) > 0, datasetName, sessionName)
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
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)

    savedPath <- paste(sessionPath, 'sess.bin', sep=.Platform$file.sep)
    if (file.exists(savedPath)) {
        return (jsonlite::toJSON(basename(savedPath)))
    }

    stop('Invalid session key')
}
