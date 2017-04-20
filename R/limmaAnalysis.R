#' @name limmaAnalysis
#' @title Differential Expression analysis
#' @description Function for performing differential expression analysis. Returns serialized de-matrix.
#' @param es an ExpressionSet object, should be normalized
#' @param columns vector of specified columns' indices (optional)
#' @param rows vector of specified rows' indices (optional)
#' @param fieldValues vector of comparison values, matching categories' names with columns (must be equal length with columns' vector if specified)
#' @return file with serialized information on differential expression matrix
#' @export
#' @import Biobase
#' @import limma
#' @import protolite
limmaAnalysis <- function(es, rows = c(), columns = c(), fieldValues) {
  assertthat::assert_that(length(columns) == length(fieldValues) || length(columns) == 0)

  rows <- getIndicesVector(rows, nrow(exprs(es)))
  columns <- getIndicesVector(columns, ncol(exprs(es)))

  fieldName <- "Comparison"
  fieldValues <- replace(fieldValues, fieldValues == '', NA)

  new.pdata <- pData(es)[columns,]
  new.pdata[[fieldName]] <- as.factor(fieldValues)
  new.pdata <- new.pdata[!is.na(new.pdata[[fieldName]]),]
  new.sampleNames <- row.names(new.pdata)

  es.copy <- es[rows, new.sampleNames]
  pData(es.copy) <- new.pdata

  fvarLabels(es.copy) <- gsub(pattern = "rowNames", replacement = "symbol", x = fvarLabels(es.copy))
  es.design <- model.matrix(~0 + pData(es.copy)[[fieldName]], data = pData(es.copy))
  colnames(es.design) <- gsub(pattern = "pData.es.copy...fieldName..",
                              replacement = '',
                              x = make.names(colnames(es.design)))

  fit <- lmFit(es.copy, es.design)
  fit2 <- contrasts.fit(fit, makeContrasts(B - A,
                                           levels=es.design))
  fit2 <- eBayes(fit2)
  de <- topTable(fit2, adjust.method="BH", number=Inf)
  de <- de[row.names(fData(es.copy)),]
  f <- tempfile(pattern = "de", tmpdir = getwd(), fileext = ".bin")
  writeBin(protolite::serialize_pb(as.list(de)), f)
  f
}

