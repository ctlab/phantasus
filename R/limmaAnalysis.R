#' Differential Expression analysis.
#'
#' \code{limmaAnalysis} performs differential expression analysis
#'   from limma package and returns a ProtoBuf-serialized resulting
#'   de-matrix.
#'
#' @param es ExpressionSet object. It should be normalized for
#'   more accurate analysis.
#'
#' @param columns Vector of specified columns' indices (optional).
#'
#' @param rows Vector of specified rows' indices (optional).
#'
#' @param fieldValues Vector of comparison values, mapping
#'   categories' names to columns/samples
#'   (must be equal length with columns' vector if specified).
#'
#' @return Name of the file containing serialized de-matrix.
#'
#' @export
#' @import Biobase
#' @import limma
limmaAnalysis <- function(es, rows = c(), columns = c(), fieldValues) {
    assertthat::assert_that(length(columns) == length(fieldValues)
                            || length(columns) == 0)

    rows <- getIndicesVector(rows, nrow(exprs(es)))
    columns <- getIndicesVector(columns, ncol(exprs(es)))

    fieldName <- "Comparison"
    fieldValues <- replace(fieldValues, fieldValues == "", NA)

    new.pdata <- pData(es)[columns, ]
    new.pdata[[fieldName]] <- as.factor(fieldValues)
    new.pdata <- new.pdata[!is.na(new.pdata[[fieldName]]), ]
    new.sampleNames <- row.names(new.pdata)

    es.copy <- es[rows, new.sampleNames]
    pData(es.copy) <- new.pdata
    fData(es.copy) <- data.frame(row.names = rownames(es.copy))

    es.design <- stats::model.matrix(~0 + Comparison, data = pData(es.copy))
    colnames(es.design) <- gsub(pattern = fieldName, replacement = "",
                                x = make.names(colnames(es.design)))

    fit <- lmFit(es.copy, es.design)
    fit2 <- contrasts.fit(fit, makeContrasts(B - A, levels = es.design))
    fit2 <- eBayes(fit2)
    de <- topTable(fit2, adjust.method = "BH", number = Inf)
    de <- de[row.names(fData(es.copy)), ]
    f <- tempfile(pattern = "de", tmpdir = getwd(), fileext = ".bin")
    writeBin(protolite::serialize_pb(as.list(de)), f)
    jsonlite::toJSON(f)
}
