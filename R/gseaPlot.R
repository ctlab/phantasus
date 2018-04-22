#' Returns path to an svg file with enrichment plot
#' @param fData list of annotation columns
#' @param fvarLabels vector of column names
#' @param rankBy name of the numeric column used for gene ranking
#' @param selectedGenes indexes of selected genes (starting from one, in the order of fData)
#' @param width width of the image (in inches)
#' @param height height of the image (in inches)
#' @return path to an svg file
#' @importFrom fgsea plotEnrichment fgsea
#' @importFrom ggplot2 ggtitle
gseaPlot <- function(fData, fVarLabels, rankBy, selectedGenes, width, height) {
    featureData <- data.frame(fData)
    colnames(featureData) <- fvarLabels
    ranks <- setNames(featureData[, rankBy], as.character(selectedGenes))

    pathway <- as.character(selectedGenes)


    fgseaRes <- fgsea(list(pathway), ranks, nperm=25000, nproc=1)

    pvalString <- sprintf("p-value %s",
    if (fgseaRes$nMoreExtreme == 0) {
        " < 1e-4"
    } else {
        sprintf(" = %g", fgseaRes$pval)
    })

    p <- plotEnrichment(pathway, ranks) + ggtitle(pvalString)
    f <- tempfile(pattern = "enrichment", tmpdir = getwd(), fileext = ".svg")
    ggsave(p, filename = f, width=width, height=height)
    jsonlite::toJSON(f)
}
