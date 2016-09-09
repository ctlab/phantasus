createES <- function(data, pData, labelDescription, colNames, rowNames) {
  stopifnot(require(Biobase))
  exprs <- t(data)
  colnames(exprs) <- colNames
  truePData <- pData
  pd <- data.frame(truePData, row.names = colNames)
  names(pd) <- labelDescription
  require(Biobase)
  phenoData <- AnnotatedDataFrame(pd)
  featureData <- AnnotatedDataFrame(data.frame(rowNames))
  es <- ExpressionSet(assayData = exprs, phenoData=phenoData, featureData = featureData)
  return(es)
}