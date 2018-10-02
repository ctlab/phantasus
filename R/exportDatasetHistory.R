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
            complexArgument <- paste(argumentName, ' <- ', toString(env[[argumentName]]))
            complexArguments <- paste(complexArguments, complexArgument, sep='\n')
        }
    }
    myHistory <- paste(complexArguments, rds[[1]]$src, sep='\n')

    if (!is.element(toString(fn), rootFn)) {
        parent <- parsed[[1]]$es[[2]]
        parentVar <- parsed[[1]]$es[[3]]
        parentHistory <- exportDatasetHistory(parent, parentVar, leaf = F)
        ourHistory <- paste(parentHistory, myHistory, sep='\n')
        if (leaf) {
            return(jsonlite::toJSON(ourHistory))
        }
        return(ourHistory)
    }
    if (leaf) {
        return(jsonlite::toJSON(myHistory))
    }
    myHistory
}
playground <- function () {
    exportDatasetHistory('x04fdf2b9ebbbda');
}
