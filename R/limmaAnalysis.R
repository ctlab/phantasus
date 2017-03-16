limmaAnalysis <- function(es, fieldName, fieldValues) {
  stopifnot(require(Biobase))
  stopifnot(require(limma))
  stopifnot(require(protolite))
  fieldValues <- replace(fieldValues, fieldValues == "", "C")
  pData(es)[[fieldName]] <- as.factor(fieldValues)

  fvarLabels(es) <- gsub(pattern = "rowNames", replacement = "symbol", x = fvarLabels(es))
  es.design <- model.matrix(~0 + pData(es)[[fieldName]], data = pData(es))
  colnames(es.design) <- gsub(pattern = "pData.es...fieldName..",
                              replacement = '',
                              x = make.names(colnames(es.design)))

  # A <- make.names(paste(field, classA, sep = ''))
  # B <- make.names(paste(field, classB, sep = ''))
  # colnames(es.design) <- gsub(pattern = A, replacement = "A", x = colnames(es.design))
  # colnames(es.design) <- gsub(pattern = B, replacement = "B", x = colnames(es.design))

  fit <- lmFit(es, es.design)
  fit2 <- contrasts.fit(fit, makeContrasts(B - A,
                                           levels=es.design))
  fit2 <- eBayes(fit2)
  de <- topTable(fit2, adjust.method="BH", number=Inf)
  de[["symbol"]] <- as.character(de[["symbol"]])
  de <- de[row.names(fData(es)),]
  f <- tempfile(pattern = "de", tmpdir = getwd(), fileext = ".bin")
  writeBin(protolite::serialize_pb(as.list(de)), f)
  f
}
