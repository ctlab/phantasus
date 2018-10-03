rootFn <- c('loadGEO', 'loadPreloaded', 'createES')
ignoredArgs <- c('es')

exportDatasetHistory <- function (sessionName, esVariable="es", leaf = T) {
    ocpuRoot <- strsplit(getwd(), 'ocpu-temp')[[1]][1]
    sessionPath <- paste(ocpuRoot, 'ocpu-store', sessionName, sep=.Platform$file.sep)
    RDataPath <- paste(sessionPath, '.RData', sep=.Platform$file.sep)
    REvalPath <- paste(sessionPath, '.REval', sep=.Platform$file.sep)
    rds <- readRDS(REvalPath)
    env <- new.env()
    loadedVariables <- load(RDataPath, envir = env)
    parsed <- parse(text=rds[[1]]$src)
    fn <- parsed[[1]][[1]]
    complexArguments <- "";
    for(i in 2:length(parsed[[1]])) {
        argumentName <- toString(parsed[[1]][[i]])
        if (argumentName %in% loadedVariables) {
            rawArgument <- paste(deparse(env[[argumentName]]), collapse = '\n')
            complexArgument <- paste(argumentName, ' <- ', rawArgument)
            complexArguments <- paste(complexArguments, complexArgument, sep='\n')
        }
    }

    if ('es' %in% names(parsed[[1]])) {
        parsedFnCall <- parsed[[1]]
        parsedFnCall$es <- parsedFnCall$es[[3]]
        rawFnCall <- paste(deparse(parsedFnCall), collapse = '')
        myHistory <- paste(complexArguments, rawFnCall, sep="\n")
    } else {
        myHistory <- paste(complexArguments, rds[[1]]$src, sep='\n')
    }

    if (!is.element(toString(fn), rootFn)) {
        parent <- parsed[[1]]$es[[2]]
        parentVar <- parsed[[1]]$es[[3]]
        parentHistory <- exportDatasetHistory(parent, parentVar, leaf = F)
        myHistory <- paste(parentHistory, myHistory, sep='\n')
    }

    if (leaf) {
        myHistory <- paste('library(phantasus)', myHistory, sep = '\n')
        return(jsonlite::toJSON(myHistory))
    }
    myHistory
}
playground <- function () {
    exportDatasetHistory('x0e04ac0ccc566d');
}
