probeDataset <- function (es, query) {
    exprsIndices <- query$exprs;
    fDataQuery <- query$fData;
    pDataQuery <- query$pData;
    response <- list()
    response[['dims']] <- dim(exprs(es))
    response[['fvarLabels']] <- colnames(fData(es))
    response[['varLabels']] <- colnames(pData(es))
    response[['probe']] <- exprs(es)[exprsIndices]
    response[['fdata']] <- list()
    response[['pdata']] <- list()
    if (!is.null(nrow(fDataQuery))) {
        for(i in seq_len(nrow(fDataQuery))) {
            row <- fDataQuery[i,]
            response[['fdata']][[row$name]] <- fData(es)[[row$name]][row$indices[[1]]]
        }
    }
    if (!is.null(nrow(pDataQuery))) {
        for(i in seq_len(nrow(pDataQuery))) {
            row <- pDataQuery[i,]
            response[['pdata']][[row$name]] <- pData(es)[[row$name]][row$indices[[1]]]
        }
    }
    jsonlite::toJSON(response)
}
