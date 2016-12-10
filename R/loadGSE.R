loadGSE <- function(name, type) {
  stopifnot(require(Biobase))
  stopifnot(require(GEOquery))
  #stopifnot(require(org.Mm.eg.db))
  #stopifnot(require(limma))
  #stopifnot(require(data.table))
  if (type == 'GSE') {
    es.loaded <- getGEO(name, AnnotGPL = T)[[1]]
    data <- as.matrix(exprs(es.loaded)); exprs <- data; colnames(data) <- NULL; row.names(data) <- NULL
    pd <- as.matrix(pData(es.loaded)); pData <- phenoData(es.loaded); colnames(pd) <- NULL; row.names(pd) <- NULL
    participants <- row.names(pData(es.loaded))
    rownames <- featureNames(es.loaded)
    fdata <- as.matrix(fData(es.loaded)[,grepl("symbol", varLabels(featureData(es.loaded)), ignore.case = T)]); colnames(fdata) <- NULL; row.names(fdata) <- NULL
    varlabels <- if (ncol(fdata) > 0) "symbol" else NULL
    res <- list(data = data, pdata = pd, participants = participants, symbol = fdata[,1], rownames = rownames, colMetaNames = colnames(pData(es.loaded)))
  }
  else {
    l <- getGEO(name)
    table <- slot(l, 'dataTable')
    data <- Table(table)
    rownames <- as.vector(data[["ID_REF"]])
    symbol <- as.vector(data[["IDENTIFIER"]])
    data <- data[,3:ncol(data)]
    participants <- colnames(data)
    columnsMeta <- Columns(table); row.names(columnsMeta) <- participants
    columnsMeta <- as.matrix(columnsMeta[,!(colnames(columnsMeta) %in% c('sample'))]); pData <- AnnotatedDataFrame(data.frame(columnsMeta))
    dm <- as.matrix(data); exprs <- dm; row.names(exprs) <- rownames; colnames(dm) <- NULL
    colMetaNames <- colnames(columnsMeta); colnames(columnsMeta) <- NULL
    res <- list(data = dm, pdata = columnsMeta, participants = participants, symbol = symbol, rownames = rownames, colMetaNames = colMetaNames)
  }
  fData <- data.frame(matrix(res$symbol, nrow(res$data), 1)); colnames(fData) <- "symbol"
  fData <- AnnotatedDataFrame(fData); featureNames(fData) <- res$rownames

  assign("es", ExpressionSet(assayData = exprs, phenoData = pData, featureData = fData), envir = parent.frame())
  f <- tempfile(pattern = "gse", tmpdir = getwd(), fileext = ".bin")
  writeBin(protolite::serialize_pb(res), f)
  f
}
