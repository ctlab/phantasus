sessionExists <- function(sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)

    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    if (file.exists(RDataPath)) {
        savedEnv <- load(RDataPath)

        return (jsonlite::toJSON(list(result="es" %in% savedEnv), auto_unbox = TRUE))
    }

    return (jsonlite::toJSON(list(result=FALSE), auto_unbox = TRUE))
}

loadSesssion <- function (sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)

    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    if (file.exists(RDataPath)) {
        savedEnv <- load(RDataPath)

        if (!("es" %in% savedEnv)) {
            return (jsonlite::toJSON(NULL))
        }

        result <- list(es=writeToList(es))
        f <- tempfile(pattern = "gse", tmpdir = getwd(), fileext = ".bin")
        writeBin(protolite::serialize_pb(result), f)
        return (jsonlite::toJSON(basename(f)))
    }

    return(jsonlite::toJSON(NULL))
}
