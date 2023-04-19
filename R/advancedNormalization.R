#'
#' @import Biobase
#' @import limma
#' @import edgeR
#'
tmmNormalization <- function(es, fieldName, logratioTrim, sumTrim) {
    if (fieldName != "None"){
        fieldValues <- es[[fieldName]]
    } else {
        fieldValues <- NULL
    }
    es.copy <- es
    count_factors <- DGEList(counts = exprs(es.copy), group = fieldValues)
    count_factors <- calcNormFactors(count_factors, logratioTrim = as.numeric(logratioTrim), sumTrim = as.numeric(sumTrim))
    exprs(es.copy) <- cpm(count_factors)

    assign("es", es.copy, envir = parent.frame())
    f <- tempfile(pattern = "norm_counts", tmpdir = getwd(), fileext = ".bin")
    data <- as.matrix(exprs(es.copy))
    colnames(data) <- NULL
    row.names(data) <- NULL
    rownames <- rownames(es.copy)

    writeBin(protolite::serialize_pb(list(data = data, colMetaNames = varLabels(es.copy))), f)
    return( jsonlite::toJSON(f))
}



voomNormalization <- function(es, fieldNames){
    if (length(fieldNames) == 0 ){
        design <- NULL
    }
    if (! "design" %in% ls()){
        fieldNames <- unlist(fieldNames)
        pdata <- pData(es)[,fieldNames]
        pdata[,fieldNames] <-  lapply(fieldNames, function(x){
            curFact <- as.factor(pdata[,x])
            curFact
        })
        design <- model.matrix(~ 0 + ., data = pdata)
        if(qr(design)$rank < ncol(design)){
            stop("Error: redundancy of model parameters appears. Try to exclude nested factors.")
        }
    }

    es.copy <- es
    voom_counts <- voom(counts = exprs(es.copy), design = design)
    exprs(es.copy) <- voom_counts$E

    assign("es", es.copy, envir = parent.frame())
    f <- tempfile(pattern = "voom_counts", tmpdir = getwd(), fileext = ".bin")
    data <- as.matrix(exprs(es.copy))
    colnames(data) <- NULL
    row.names(data) <- NULL
    rownames <- rownames(es.copy)

    writeBin(protolite::serialize_pb(list(data = data, colMetaNames = varLabels(es.copy))), f)
    return( jsonlite::toJSON(f))
}
