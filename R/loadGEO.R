#' @name loadGEO
#' @title Load GEO Dataset
#' @description Function for loading dataset from GEO, creating ExpressionSet from it and exporting it to global environment
#' @param name string, containing GEO identifier of the dataset
#' @param type "GSE" or "GDS"
#' @return file with serialized information on ExpressionSet
#' @export
#' @import Biobase
#' @import GEOquery
loadGEO <- function(name, type = NA) {
  assign("es", getES(name, type, destdir = "/var/morpheus/cache"), envir = parent.frame())
  data <- as.matrix(exprs(es)); colnames(data) <- NULL; row.names(data) <- NULL

  pdata <- as.matrix(pData(es)); colnames(pdata) <- NULL; row.names(pdata) <- NULL

  participants <- colnames(es)
  rownames <- rownames(es)

  fdata <- as.matrix(fData(es)[,grepl("symbol", varLabels(featureData(es)), ignore.case = T)])
  colnames(fdata) <- NULL
  row.names(fdata) <- NULL

  varlabels <- if (ncol(fdata) > 0) "symbol" else NULL
  fdata <- if (is.null(varlabels)) NULL else fdata[,1]

  res <- list(data = data, pdata = pdata, participants = participants, symbol = fdata, rownames = rownames, colMetaNames = colnames(pData(es.loaded)))

  f <- tempfile(pattern = "gse", tmpdir = getwd(), fileext = ".bin")
  writeBin(protolite::serialize_pb(res), f)
  f
}
