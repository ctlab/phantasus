rootFn <- c('loadGEO', 'loadPreloaded', 'createES')

exportDatasetHistory <- function (sessionName, esVariable="es") {
    ocpuRoot <- Sys.getenv('OCPU_MASTER_HOME')
    if (length (ocpuRoot) <= 1) {
        ocpuRoot <- tempdir()
    }
    REvalPath <- paste(ocpuRoot, 'ocpu-store',sessionName, '.REval', sep=.Platform$file.sep)
    rds <- readRDS(REvalPath)
    parsed <- parse(text=rds[[1]]$src)
    fn <- parsed[[1]][[1]]
    if (!is.element(toString(fn), rootFn)) {
        parent <- parsed[[1]]$es[[2]]
        parentVar <- parsed[[1]]$es[[3]]
        parentHistory <- exportDatasetHistory(parent, parentVar)
        ourHistory <- paste(parentHistory, rds[[1]]$src, sep='\n')
        return(jsonlite::toJSON(ourHistory))
    }
    rds[[1]]$src

}
playground <- function () {
    exportDatasetHistory('x093ddc6795');
}
