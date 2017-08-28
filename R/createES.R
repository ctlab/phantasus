#' Create ExpressionSet.
#'
#' \code{createES} function produces an ExpressionSet object from given data,
#'     and exports it to global scope.
#'
#' @param data Gene expression matrix.
#'
#' @param pData Matrix with phenotypical data.
#'
#' @param varLabels Names of phenoData columns.
#'
#' @param fData Matrix with feature data.
#'
#' @param fvarLabels Names of featureData columns.
#'
#' @return produced ExpressionSet object
#'
#' @export
#' @import Biobase
#'
#' @examples
#' data <- matrix(1:15, 5, 3)
#' pData <- c("A", "B", "C")
#' varLabels <- "cat"
#' fData <- c("p", "r", "s", "t", "u")
#' fvarLabels <- "id"
#' createES(data, pData, varLabels, fData, fvarLabels)
#'
createES <- function(data, pData, varLabels, fData, fvarLabels) {
    phenoData <- AnnotatedDataFrame(data.frame(pData))
    varLabels(phenoData) <- varLabels

    featureData <- AnnotatedDataFrame(data.frame(fData))
    varLabels(featureData) <- fvarLabels

    es <- ExpressionSet(assayData = data,
                        phenoData = phenoData,
                        featureData = featureData)
    assign("es", es, envir = parent.frame())
    es
}
