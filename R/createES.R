#' @name createES
#' @title Create ExpressionSet
#' @description Function for creating ExpressionSet from parameters and exporting it to global environment
#' @param data gene expression matrix
#' @param pData matrix with phenotypical data
#' @param labelDescription names of phenoData
#' @param colNames sample ids
#' @param rowNames gene ids
#' @param symbol gene symbol
#' @export
#' @import Biobase
createES <- function(data, pData, varLabels, fData, fvarLabels) {
  exprs <- t(data)
  colnames(exprs) <- colNames
  truePData <- pData
  phenoData <- AnnotatedDataFrame(data.frame(pData))
  varLabels(phenoData) <- varLabels
  
  featureData <- AnnotatedDataFrame(data.frame(fData))
  varLabels(featureData) <- fvarLabels
 
  es <- ExpressionSet(assayData = exprs, phenoData=phenoData, featureData = featureData)
  assign("es", es, envir = parent.frame())
  es
}
