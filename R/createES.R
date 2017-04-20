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
createES <- function(data, pData, labelDescription, colNames, rowNames, symbol) {
  exprs <- t(data)
  colnames(exprs) <- colNames
  truePData <- pData
  pd <- data.frame(truePData, row.names = colNames)
  names(pd) <- labelDescription
  phenoData <- AnnotatedDataFrame(pd)
  if (is.null(symbol) || length(symbol) == 0) {
    symbol <- 1:nrow(exprs)
  }
  featureData <- AnnotatedDataFrame(data.frame(symbol))
  if (is.null(rowNames) || length(rowNames) == 0) {
    rowNames <- 1:nrow(exprs)
  }
  featureNames(featureData) <- rowNames
  es <- ExpressionSet(assayData = exprs, phenoData=phenoData, featureData = featureData)
  assign("es", es, envir = parent.frame())
  es
}
