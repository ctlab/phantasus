#' Collapse dataset
#'
#' \code{collapseDataset} performs a collapse action on expression set
#'
#' @param es Expression set
#' @param isRows Work with rows. False if columns (default True - row mode)
#' @param selectOne select best match or merge duplicates
#' @param fn select/merge function
#' @param fields fields to unique on
#'
#' @export
#' @import ccaPP
#'
#' @examples
#' \dontrun{
#' es <- getGSE('GSE53986')[[1]]
#' collapseDataset(es, isRows = TRUE, selectOne = TRUE,
#' fn = mean, fields = c('Gene ID', 'Gene symbol'))
#' }
#'
collapseDataset <- function (es, isRows = TRUE, selectOne = FALSE, fn, fields) {
    es <- collapseDatasetImpl(es, isRows, selectOne, fn, fields)

    assign("es", es, envir = parent.frame())
}

collapseDatasetImpl <- function (es, isRows = TRUE, selectOne = FALSE, fn, fields) {
    expr <- exprs(es)
    fact <- collectFactor(es, isRows, fields)
    f2 <- factor(fact, levels=unique(fact))

    if (selectOne) { #select fittest
        ranks <- apply(expr, 1, fn)
        factorFrame <- data.frame(f=f2, i=seq_along(ranks), r=ranks)
        factorFrame <- factorFrame[order(factorFrame$f, -factorFrame$r), ]
        keep <- factorFrame[!duplicated(factorFrame$f) & !is.na(factorFrame$f), ]$i
        res <- es[keep, ]
        return(res)
    } else { #merge duplicates
        factorFrame <- data.frame(f=f2, i=seq_len(length(f2)))
        keep <- factorFrame[!duplicated(factorFrame$f) & !is.na(factorFrame$f), ]$i
        keep <- sort(keep)

        if (isRows) {
            res <- es[keep, ]
            oldAnnotation <- fData(res)
        } else {
            expr <- t(expr)
            res <- es[, keep]
            oldAnnotation <- pData(res)
        }

        splitted <- split(seq_len(length(f2)), fact)
        zz <- lapply(seq_len(ncol(expr)), function(i) split(expr[, i], f2))
        zz <- unlist(zz, recursive = FALSE)
        collapsedExprs <- lapply(seq_len(ncol(expr)),
                            function (i) sapply(split(expr[, i], f2), fn))
        collapsedExprs <- do.call(cbind, collapsedExprs)

        if (!isRows) {
            collapsedExprs <- t(collapsedExprs)
        }

        rownames(collapsedExprs) <- rownames(res)
        colnames(collapsedExprs) <- colnames(res)
        exprs(res) <- collapsedExprs
        fields <- colnames(oldAnnotation)[which(colnames(oldAnnotation) %in% fields)]
        newAnnotaion <- oldAnnotation[, which(colnames(oldAnnotation) %in% fields), drop=F]
        rownames(newAnnotaion) <- rownames(oldAnnotation)
        colnames(newAnnotaion) <- fields
        if (isRows) {
            fData(res) <- newAnnotaion
        } else {
            pData(res) <- newAnnotaion
        }

        return(res)
    }
}

collectFactor <- function (es, isRows, fields) {
    if (!length(fields)) {
        stop('Empty fields given')
    }

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
