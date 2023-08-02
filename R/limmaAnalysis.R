
limmaAnalysisSimpleImpl <- function(es, fieldValues, contrast){
    fieldValues <- replace(fieldValues, fieldValues == "", NA)

    es.copy <- es
    es.copy$Comparison <- fieldValues
    pData(es.copy)[,"Comparison"] <- as.factor(pData(es.copy)[,"Comparison"])
    pData(es.copy)[,"Comparison"] <- relevel(pData(es.copy)[,"Comparison"], ref = "Reference")
    fData(es.copy) <- data.frame(row.names = rownames(es.copy))
    es.copy <- es.copy[, !is.na(fieldValues)]

    # Getting rid of check NOTEs
    Comparison=ComparisonA=ComparisonB=NULL

    es.design <- stats::model.matrix(~0 + Comparison, data = pData(es.copy))

    fit <- lmFit(es.copy, es.design)

    A <- NULL; B <- NULL
    fit2 <- contrasts.fit(fit, makeContrasts(ComparisonTarget - ComparisonReference, levels = es.design))
    fit2 <- eBayes(fit2)
    de <- topTable(fit2, adjust.method = "BH", number = Inf)
    de <- de[row.names(fData(es.copy)), ]
    return(de)
}
limmaAnalysisAdvancedImpl <- function(es, designData, contrast){
    ux_designMatrix <- getDesignMatrix(designData)

    es.copy <- es
    fData(es.copy) <- data.frame(row.names = rownames(es.copy))
    colnames(ux_designMatrix) <-  make.names(colnames(ux_designMatrix))
    target_level <- make.names(paste0(contrast[1],contrast[3]))
    reference_level <-  make.names(paste0(contrast[1],contrast[2]))


    fit <- lmFit(es.copy, ux_designMatrix)
    fit2 <- contrasts.fit(fit, makeContrasts(contrasts = paste(reference_level,target_level, sep = "-"), levels = ux_designMatrix))
    fit2 <- eBayes(fit2)
    de <- topTable(fit2, adjust.method = "BH", number = Inf)
    de <- de[row.names(fData(es.copy)), ]
    return(de)
}


#' Differential Expression analysis.
#'
#' \code{limmaAnalysis} performs differential expression analysis
#'     from limma package and returns a ProtoBuf-serialized resulting
#'     de-matrix.
#'
#' @param es ExpressionSet object. It should be normalized for
#'     more accurate analysis.
#'
#' @param fieldValues Vector of comparison values, mapping
#'     categories' names to columns/samples
#'
#' @param version name of the limma analysis implementation. Should be "One-factor design" or "Advanced design"
#'
#' @param contrast a character vector with exactly three elements: the name of a factor in the design formula, the name of the numerator level for the fold change, and the name of the denominator level for the fold change
#'
#' @param designData data.frame with design matrix
#'
#' @return Name of the file containing serialized de-matrix.
#'
#' @import Biobase
#' @import limma
#'
#' @examples
#' \dontrun{
#' data(es)
#' limmaAnalysis(es, fieldValues = c("A", "A", "A", "B", "B"))
#' }
limmaAnalysis <- function (es, fieldValues, version = "One-factor design", contrast =  list('Comparison', 'Target', 'Reference'), designData = NULL) {
    fieldValues <- replace(fieldValues, fieldValues == "", NA)
    de <- NULL
    contrast <- unlist(contrast)
    if (version == "One-factor design" ){
        de <- limmaAnalysisSimpleImpl(es, fieldValues, contrast)
    }
    if (version == "Advanced design"){
        de <- limmaAnalysisAdvancedImpl(es, designData, contrast)
    }
        deDf <- as.data.frame(de)
    toRemove <- intersect(colnames(fData(es)), colnames(deDf))
    fData(es)[, toRemove] <- NULL

    es$Comparison <- fieldValues
    fData(es) <- cbind(fData(es), deDf)
    assign("es", es, envir = parent.frame())

    f <- tempfile(pattern = "de", tmpdir = getwd(), fileext = ".bin")
    writeBin(protolite::serialize_pb(as.list(de)), f)
    jsonlite::toJSON(f)
}

