#' @name createES
#' @title Create ExpressionSet
#' @description Function for creating ExpressionSet from parameters and exporting it to global environment
#' @param data gene expression matrix
#' @param pData matrix with phenotypical data
#' @param varLabels names of phenoData columns
#' @param fData matrix with feature data
#' @param fvarLabels names of featureData columns
#' @export
#' @import Biobase
createES <- function(data, pData, varLabels, fData, fvarLabels) {
    phenoData <- AnnotatedDataFrame(data.frame(pData))
    varLabels(phenoData) <- varLabels

    featureData <- AnnotatedDataFrame(data.frame(fData))
    varLabels(featureData) <- fvarLabels

    es <- ExpressionSet(assayData = data,
                        phenoData = phenoData,
                        featureData = featureData)
    assign("es", es, envir = parent.frame())
    es
}
