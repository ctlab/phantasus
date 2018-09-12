collapseDatasetImpl <- function (es, isRows = TRUE, selectOne = FALSE, fn, fields) {
    expr <- exprs(es)
    fact <- collectFactor(es, isRows, fields)
    f2 <- factor(fact, levels=unique(fact))

    if (!selectOne) {
        t <- data.frame(f=f2, i=seq_len(length(f2)))
        keep <- t[!duplicated(t$f) & !is.na(t$f),]$i
        keep <- sort(keep)
        if (isRows) {
            res <- es[keep, ]
            splitted <- split(seq_len(length(f2)), fact)
            zz <- lapply(seq_len(ncol(expr)), function(i) split(expr[,i], f2))
            zz <- unlist(zz, recursive = FALSE)
            collapsed <- lapply(seq_len(ncol(expr)),
                                function (i) sapply(split(expr[,i], f2), fn))
            collapsed <- do.call(cbind, collapsed)
            rownames(collapsed) <- rownames(res)
            colnames(collapsed) <- colnames(res)
            exprs(res) <- collapsed
            newFdata <- as.data.frame(fData(res)[, which(colnames(fData(res)) %in% fields)])
            rownames(newFdata) <- rownames(fData(res))
            colnames(newFdata) <- fields
            fData(res) <- newFdata
            return(res)
        } else {
            #not working
            res <- es[, keep]
            splitted <- split(seq_len(length(f2)), fact)
            zz <- lapply(seq_len(nrow(expr)), function(i) split(expr[i,], f2))
            zz <- unlist(zz, recursive = FALSE)
            collapsed <- lapply(seq_len(nrow(expr)),
                                function (i) sapply(split(expr[i,], f2), fn))
            collapsed <- do.call(cbind, collapsed)
            collapsed <- t(collapsed)
            rownames(collapsed) <- rownames(res)
            colnames(collapsed) <- colnames(res)
            exprs(res) <- collapsed
            newPhenoData <- as.data.frame(pData(res)[, which(colnames(pData(res)) %in% fields)])
            rownames(newPhenoData) <- rownames(phenoData(res))
            colnames(newPhenoData) <- fields
            pData(res) <- newPhenoData
            return(res)
        }
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

collectFactor <- function (es, isRows, fields) {
    if (isRows) {
        target <- fData(es)
    } else {
        target <- phenoData(es)
    }


    f <- target[[fields[1]]]
    if (length(fields) > 1) {
        for(i in 2:length(fields)) {
            f <- paste(f, target[[fields[i]]], sep='//r')
        }
    }
    return(f)
}

collapseDataset <- function (es, isRows = TRUE, selectOne = FALSE, fn, fields) {
    es <- collapseDatasetImpl(es, isRows, selectOne, fn, fields)

    assign("es", es, envir = parent.frame())
}

playground <- function () {
    es <- getGSE("GSE53986")[[1]]
    collapseDataset(es,selectOne = FALSE, fn = median, fields = c('Gene ID', 'Gene symbol'))
    collapseDataset(es, isRows = FALSE, selectOne = FALSE, fn = median, fields = c('treatment'))
}
