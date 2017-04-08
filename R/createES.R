createES <- function(data, pData, labelDescription, colNames, rowNames, symbol) {
  stopifnot(require(Biobase))
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
  return(es)
}
