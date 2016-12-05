loadGSE <- function(name) {
  stopifnot(require(Biobase))
  stopifnot(require(GEOquery))
  #stopifnot(require(org.Mm.eg.db))
  #stopifnot(require(limma))
  #stopifnot(require(data.table))
  es <- getGEO(name)[[1]]
  res <- list(data = as.vector(exprs(es)), dimData = c(ncol(exprs(es)), nrow(exprs(es))), pData = as.vector(as.matrix(pData(es))),
              dimPData = c(ncol(pData(es)), nrow(pData(es))),)



#
#   pData(es)$condition <- sub("-.*$", "", es$title)
#
#   es <- collapseBy(es, fData(es)$ENTREZ_GENE_ID)
#   es <- es[!grepl("///", rownames(es)), ]
#   es <- es[rownames(es) != "", ]
#
#   # there is a lot of garbage there
#   fData(es) <- data.frame(row.names = rownames(es))
#   fData(es)$symbol <- sapply(mget(rownames(es),
#                                   org.Mm.egSYMBOL,
#                                   ifnotfound = NA),
#                              head, n=1)
#   exprs(es) <- normalizeBetweenArrays(log2(exprs(es) + 1), method="quantile")
#   es.design <- model.matrix(~0+condition, data=pData(es))
  return(protolite::serialize_pb(es))
}

collapseBy <- function(es, factor, FUN=median) {
  ranks <- apply(exprs(es), 1, FUN)
  t <- data.frame(f=factor, i=seq_along(ranks), r=ranks)
  t <- t[order(t$r, decreasing=T), ]
  keep <- t[!duplicated(t$f) & !is.na(t$f),]$i
  res <- es[keep, ]
  fData(res)$origin <- rownames(res)
  rownames(res) <- factor[keep]
  res
}
