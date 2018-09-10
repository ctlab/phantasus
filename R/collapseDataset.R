collapseDatasetImpl <- function (es, isRows = TRUE, selectOne = FALSE, fn, fields) {
    expr <- exprs(es)
    f <- fData(es)[[fields[1]]]
    for(i in 2:length(fields)) {
        f <- paste(f, fData(es)[[fields[i]]], sep='//r')
    }
    f2 <- factor(f, levels=unique(f))

    if (!selectOne) {
        t <- data.frame(f=f2, i=seq_along(expr))
        keep <- t[!duplicated(t$f) & !is.na(t$f),]$i
        keep <- sort(keep)
        res <- es[keep, ]
        #rownames(res) <- f2[keep]


        splitted <- split(seq_len(length(f2)), f)
        zz <- lapply(seq_len(ncol(expr)), function(i) split(expr[,i], f2))
        zz <- unlist(zz, recursive = FALSE)
        collapsed <- lapply(seq_len(ncol(expr)),
                            function (i) sapply(split(expr[,i], f2), fn))
        collapsed <- do.call(cbind, collapsed)
        rownames(collapsed) <- rownames(res)
        exprs(res) <- collapsed
        if (isRows) {
            fData(res) <- fData(res)[, which(colnames(fData(res)) %in% fields)]
        } else {
            phenoData(res) <- phenoData(res)[, which(colnames(phenoData(res)) %in% fields)]
        }
        browser()
    } else {
        ranks <- apply(expr, 1, fn)
        t <- data.frame(f=f2, i=seq_along(ranks), r=ranks)
        t <- t[order(t$r, decreasing=T), ]
        keep <- t[!duplicated(t$f) & !is.na(t$f),]$i
        keep <- sort(keep)
        res <- es[keep, ]
        #rownames(res) <- f2[keep]
        return(res)
    }
}

collapseDataset <- function (es, isRows = TRUE, selectOne = FALSE, fn, fields) {
    es <- collapseDatasetImpl(es, isRows, selectOne, fn, fields)

    assign("es", es, envir = parent.frame())
}

playground <- function () {
    es <- getGSE("GSE53986")[[1]]
    collapseDataset(es,selectOne = FALSE, fn = median, fields = c('Gene ID', 'Gene symbol'))
}
