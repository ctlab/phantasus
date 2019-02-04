sessionExists <- function(sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)

    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    if (file.exists(RDataPath)) {
        env <- new.env()
        load(RDataPath, envir = env)
        validForLoading = !is.null(env$es) && !is.null(attr(env$es, 'published'))

        return (jsonlite::toJSON(list(result=validForLoading), auto_unbox = TRUE))
    }

    return (jsonlite::toJSON(list(result=FALSE), auto_unbox = TRUE))
}

publishSession <- function (sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)

    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    if (file.exists(RDataPath)) {
        env <- new.env()
        load(RDataPath, envir = env)

        if (!is.null(env$es)) { # this session is for actual dataset
            attr(env$es, 'published') <- TRUE
            save(list = ls(all.names = TRUE, env), file = file.path(sessionPath, ".RData"), envir = env)
            return (jsonlite::toJSON(list(result=TRUE), auto_unbox = TRUE))
        }
    }

    return (jsonlite::toJSON(list(result=FALSE), auto_unbox = TRUE))
}

loadSession <- function (sessionName) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)

    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    if (file.exists(RDataPath)) {
        savedEnv <- load(RDataPath)

        if (!("es" %in% savedEnv) || is.null(attr(es, 'published'))) {
            stop('Invalid session key')
        }

        result <- list(es=writeToList(es))
        f <- tempfile(pattern = "session", tmpdir = getwd(), fileext = ".bin")
        writeBin(protolite::serialize_pb(result), f)
        return (jsonlite::toJSON(basename(f)))
    }

    stop('Invalid session key')
}
