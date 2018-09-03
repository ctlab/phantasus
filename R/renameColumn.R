renameColumn <- function (es, isFeature = TRUE, oldName, newName) {
    if (isFeature) {
        fData(es)[[newName]] <- fData(es)[[oldName]]
        fData(es)[[oldName]] <- NULL
    } else {
        phenoData(es)[[newName]] <- phenoData(es)[[oldName]]
        phenoData(es)[[oldName]] <- NULL
    }

    assign("es", es, envir = parent.frame())
}
