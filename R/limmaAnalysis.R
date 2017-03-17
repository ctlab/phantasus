limmaAnalysis <- function(es, fieldName, fieldValues) {
  stopifnot(require(Biobase))
  stopifnot(require(limma))
  stopifnot(require(protolite))
  fieldValues <- replace(fieldValues, fieldValues == "", NA)

  new.pdata <- pData(es)
  new.pdata[[fieldName]] <- as.factor(fieldValues)
  new.pdata <- new.pdata[!is.na(new.pdata[[fieldName]]),]
  new.sampleNames <- row.names(new.pdata)
  new.exprs <- exprs(es)[,new.sampleNames]

  es.copy <- ExpressionSet(new.exprs, AnnotatedDataFrame(new.pdata), featureData = featureData(es))

  fvarLabels(es.copy) <- gsub(pattern = "rowNames", replacement = "symbol", x = fvarLabels(es.copy))
  es.design <- model.matrix(~0 + pData(es.copy)[[fieldName]], data = pData(es.copy))
  colnames(es.design) <- gsub(pattern = "pData.es.copy...fieldName..",
                              replacement = '',
                              x = make.names(colnames(es.design)))

  # A <- make.names(paste(field, classA, sep = ''))
  # B <- make.names(paste(field, classB, sep = ''))
  # colnames(es.design) <- gsub(pattern = A, replacement = "A", x = colnames(es.design))
  # colnames(es.design) <- gsub(pattern = B, replacement = "B", x = colnames(es.design))

  fit <- lmFit(es.copy, es.design)
  fit2 <- contrasts.fit(fit, makeContrasts(B - A,
                                           levels=es.design))
  fit2 <- eBayes(fit2)
  de <- topTable(fit2, adjust.method="BH", number=Inf)
  de[["symbol"]] <- as.character(de[["symbol"]])
  de <- de[row.names(fData(es.copy)),]
  f <- tempfile(pattern = "de", tmpdir = getwd(), fileext = ".bin")
  writeBin(protolite::serialize_pb(as.list(de)), f)
  f
}
