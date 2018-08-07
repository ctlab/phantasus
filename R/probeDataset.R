probeDataset <- function (es, indices) {
    response <- list()
    response[['dims']] <- dim(exprs(es))
    response[['fvarLabels']] <- colnames(fData(es))
    response[['probe']] <- exprs(es)[indices]

    jsonlite::toJSON(response)
}
