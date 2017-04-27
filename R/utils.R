getIndicesVector <- function(current, neededLength) {
  if (length(current) == 0) {
    current <- 0:(neededLength - 1)
  }
  current + 1
}

replacenas <- function(data, replacena) {
  ind <- which(is.na(data), arr.ind = T)
  if (nrow(ind) > 0) {
    data[ind] <- apply(data, 1, replacena, na.rm = T)[ind[,1]]
  }
  ind1 <- which(!is.nan(as.matrix(data)), arr.ind = T)
  left.rows <- unique(ind1[,"row"])
  data <- data[left.rows,]
  data
}

getGDS <- function(name, destdir = tempdir()) {
  l <- getGEO(name, destdir = destdir)
  table <- slot(l, 'dataTable') # extracting all useful information on dataset
  data <- Table(table)  # extracting table ID_REF | IDENTIFIER/SAMPLE | SAMPLE1 | ...
  columnsMeta <- Columns(table) # phenoData
  sampleNames <- as.vector(columnsMeta[["sample"]])
  rownames <- as.vector(data[["ID_REF"]])
  symbol <- as.vector(data[["IDENTIFIER"]])

  data <- data[,sampleNames] # expression data
  exprs <- as.matrix(data)
  row.names(exprs) <- rownames

  row.names(columnsMeta) <- sampleNames
  columnsMeta <- columnsMeta[,!(colnames(columnsMeta) %in% c('sample'))] 
  pData <- AnnotatedDataFrame(data.frame(columnsMeta, check.names = F))

  fData <- data.frame(matrix(symbol, nrow(exprs), 1));
  colnames(fData) <- "symbol"
  fData <- AnnotatedDataFrame(fData)
  featureNames(fData) <- res$rownames

  ExpressionSet(assayData = exprs, phenoData = pData, featureData = fData)
}

getES <- function(name, type = NA, destdir = tempdir()) {
  if (is.na(type)) {
     type = substr(name, 1, 3)
  }

  if (type == 'GSE') {
    es <- getGEO(name, AnnotGPL = T, destdir = destdir)[[1]]
  }
  else if (type == "GDS") {
    es <- getGDS(name, destdir)
  }
  else {
    stop("Incorrect name or type of the dataset")
  }
  es
}