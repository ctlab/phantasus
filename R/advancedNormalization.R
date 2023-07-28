#'
#' @import Biobase
#' @import limma
#' @import edgeR
#'
tmmNormalization <- function(es, fieldName, logratioTrim, sumTrim, convertCPM = FALSE) {
    if (fieldName != "None"){
        fieldValues <- es[[fieldName]]
    } else {
        fieldValues <- NULL
    }
    es.copy <- es
    count_factors <- DGEList(counts = exprs(es.copy), group = fieldValues)
    count_factors <- calcNormFactors(count_factors, logratioTrim = as.numeric(logratioTrim), sumTrim = as.numeric(sumTrim))
    if (convertCPM){
        count_factors <- cpm(count_factors)
    }
    exprs(es.copy) <- count_factors

    assign("es", es.copy, envir = parent.frame())
    f <- tempfile(pattern = "norm_counts", tmpdir = getwd(), fileext = ".bin")
    data <- as.matrix(exprs(es.copy))
    colnames(data) <- NULL
    row.names(data) <- NULL
    rownames <- rownames(es.copy)

    writeBin(protolite::serialize_pb(list(data = data, colMetaNames = varLabels(es.copy))), f)
    return( jsonlite::toJSON(f))
}


#'
#' @import Biobase
#' @import limma
#' @import edgeR
#'
voomNormalization <- function(es, designData, filterByExp = FALSE){
    designMatrix <- getDesignMatrix(designData)
    es.copy <- es
    keep <- TRUE
    if (filterByExp){
        keep <- filterByExpr(exprs(es.copy), designMatrix)
    }
    es.copy <- es.copy[keep,]
    voom_counts <- voom(counts = exprs(es.copy), design = designMatrix)

    exprs(es.copy) <- voom_counts$E

    assign("es", es.copy, envir = parent.frame())
    f <- tempfile(pattern = "voom_counts", tmpdir = getwd(), fileext = ".bin")
    data <- as.matrix(exprs(es.copy))
    colnames(data) <- NULL
    row.names(data) <- NULL

    writeBin(protolite::serialize_pb(list(data = data, keep = which(keep)-1)), f)
    return( jsonlite::toJSON(f))
}
