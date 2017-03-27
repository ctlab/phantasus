createES <- function(data, pData, labelDescription, colNames, rowNames, symbol) {
  stopifnot(require(Biobase))
  exprs <- t(data)
  colnames(exprs) <- colNames
  truePData <- pData
  pd <- data.frame(truePData, row.names = colNames)
  names(pd) <- labelDescription
  phenoData <- AnnotatedDataFrame(pd)
  featureData <- AnnotatedDataFrame(data.frame(symbol))
  featureNames(featureData) <- rowNames
  es <- ExpressionSet(assayData = exprs, phenoData=phenoData, featureData = featureData)
  assign("es", es, envir = parent.frame())
  return(es)
}
