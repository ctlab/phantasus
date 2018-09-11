probeDataset <- function (es, query) {
    exprsIndices <- query$exprs;
    fDataQuery <- query$fData;

    response <- list()
    response[['dims']] <- dim(exprs(es))
    response[['fvarLabels']] <- colnames(fData(es))
    response[['probe']] <- exprs(es)[exprsIndices]
    response[['fdata']] <- list()
    for(i in 1:nrow(fDataQuery)) {
        row <- fDataQuery[i,]
        response[['fdata']][[row$name]] <- fData(es)[[row$name]][row$indices[[1]]]
    }

    jsonlite::toJSON(response)
}
