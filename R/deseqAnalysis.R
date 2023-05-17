#https://support.bioconductor.org/p/67600/#67612
#https://github.com/ctlab/sysbio-training/blob/master/masters-2019/rnaseq/mnt/scripts/do_deseq2.R


deseqAnalysis <- function (es, fieldValues, version = "One-factor design", contrast =  list('Comparison', 'A', 'B'), designFields = NULL, designData = NULL) {
    save(es,fieldValues, version, contrast, designFields, designData, file =  "~/test_annot/deseq_param.rda")
    fieldValues <- replace(fieldValues, fieldValues == "", NA)
    de <- NULL
    contrast <- unlist(contrast)
    if (version == "One-factor design" ){
        de <- deseqAnalysisSimpleImpl(es, fieldValues, contrast)
    }
    if (version == "Advanced design"){
        de <- deseqAnalysisAdvancedImpl(es, designData, designFields, contrast)
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

deseqAnalysisSimpleImpl <- function(es, fieldValues, contrast){
    fieldValues <- replace(fieldValues, fieldValues == "", NA)

    es.copy <- es
    es.copy$Comparison <- fieldValues
    pData(es.copy)[,"Comparison"] <- as.factor(pData(es.copy)[,"Comparison"])
    pData(es.copy)[,"Comparison"] <- relevel(pData(es.copy)[,"Comparison"], ref = "B")
    es.copy <- es.copy[, !is.na(fieldValues)]
    dds <- DESeq2::DESeqDataSetFromMatrix(exprs(es.copy), pData(es.copy), design=~Comparison)

    populatedDds <- DESeq2::DESeq(dds)
    de <- DESeq2::results(populatedDds, contrast = contrast, cooksCutoff = FALSE)
    shr_values <- DESeq2::lfcShrink(populatedDds, coef =  paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_"), type = "apeglm")
    de$log2FoldChange <- shr_values$log2FoldChange
    de$lfcSE <- shr_values$lfcSE
    return(de)
}
deseqAnalysisAdvancedImpl <- function(es, designData, designFields, contrast){
    ux_designMatrix <- getDesignMatrix(designData)

    es.copy <- es
    designFields <- unlist(designFields)
    designFields <- c("~1", designFields)
    designFormula <- formula(paste(designFields, collapse = "+"))

    pData(es.copy)[,contrast[1]] <- as.factor(pData(es.copy)[,contrast[1]])
    pData(es.copy)[,contrast[1]] <- relevel(pData(es.copy)[,contrast[1]], ref = contrast[3])
    dds <- DESeq2::DESeqDataSetFromMatrix(exprs(es.copy), pData(es.copy), design=designFormula)
    real_designMatrix <- model.matrix(dds@design, data = pData(es.copy))
    real_designMatrix <- real_designMatrix[, colnames(ux_designMatrix)]
    if (!all.equal(ux_designMatrix[,], real_designMatrix[,])){
        stop("Internal error! Clientside and serverside designs are not equal!")
    }
    populatedDds <- DESeq2::DESeq(dds)
    de <- DESeq2::results(populatedDds, contrast = contrast, cooksCutoff = FALSE)
    coef_name <- paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_")
    coef_name <- gsub(pattern = "[^[:alnum:]_.]", replacement = ".", x = coef_name )
    shr_values <- DESeq2::lfcShrink(populatedDds, coef =coef_name, type = "apeglm")
    de$log2FoldChange <- shr_values$log2FoldChange
    de$lfcSE <- shr_values$lfcSE
    return(de)
}

