rootFn <- c('loadGEO', 'loadPreloaded', 'createES')
complexArgumentLimitBytes <- 400

#' Reproduce session in R code
#'
#' @param sessionName String, OCPU session name
#' @param leaf Boolean, is it leaf (default = F)
#' @param step Integer, step of recursion (default = 0)
#' @param savedEnv Environment, where to store complex arguments (default = new.env())
#'
#' @return JSON with R code
#' @importFrom utils object.size
#' @examples
#' \dontrun{
#'   setwd(tempdir())
#'   reproduceInR('x039f1672026678');
#' }
reproduceInR <- function (sessionName, leaf = T, step = 0, savedEnv = new.env()) {
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
            argumentValue <- env[[argumentName]]
            if (object.size(argumentValue) > complexArgumentLimitBytes) {
                rawArgument <- paste(argumentName, '_', step, sep='')
                savedEnv[[rawArgument]] <- argumentValue
            } else {
                rawArgument <- paste(deparse(argumentValue), collapse = '\n')
            }

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
        parentHistory <- reproduceInR(parent, leaf = F, step = step + 1, savedEnv = savedEnv)
        myHistory <- paste(parentHistory, myHistory, sep='\n')
    }

    if (leaf) {
        assign('env', savedEnv, envir = parent.frame())
        myHistory <- paste('library(phantasus)', 'load(\'~/Downloads/env.rda\') # Please note', myHistory, sep = '\n')
        return(jsonlite::toJSON(myHistory))
    }
    myHistory
}
