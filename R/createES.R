#' Create ExpressionSet.
#'
#' \code{createES} function produces an ExpressionSet object from given data,
#'   and exports it to global scope.
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
