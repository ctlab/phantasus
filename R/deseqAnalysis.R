#https://support.bioconductor.org/p/67600/#67612
#https://github.com/ctlab/sysbio-training/blob/master/masters-2019/rnaseq/mnt/scripts/do_deseq2.R

deseqAnalysis <- function (es, fieldValues) {
    fieldValues <- replace(fieldValues, fieldValues == "", NA)

    es.copy <- es
    es.copy$Comparison <- fieldValues
    es.copy <- es.copy[, !is.na(fieldValues)]
    dds <- DESeq2::DESeqDataSetFromMatrix(exprs(es.copy), pData(es.copy), design=~Comparison)

    populatedDds <- DESeq2::DESeq(dds)
    de <- DESeq2::lfcShrink(populatedDds, contrast = c('Comparison', 'A', 'B'), cooksCutoff = FALSE)

    deDf <- as.data.frame(de)
    toRemove <- intersect(colnames(fData(es)), colnames(deDf))
    fData(es)[, toRemove] <- NULL

    es$Comparsion <- fieldValues
    fData(es) <- cbind(fData(es), deDf)
    assign("es", es, envir = parent.frame())

    f <- tempfile(pattern = "de", tmpdir = getwd(), fileext = ".bin")
    writeBin(protolite::serialize_pb(as.list(de)), f)
    jsonlite::toJSON(f)
}
