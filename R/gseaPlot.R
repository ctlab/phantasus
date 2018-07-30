#' Returns path to an svg file with enrichment plot
#' @param es ExpressionSet object.
#' @param rankBy name of the numeric column used for gene ranking
#' @param selectedGenes indexes of selected genes (starting from one, in the order of fData)
#' @param width width of the image (in inches)
#' @param height height of the image (in inches)
#' @param vertical whether to use vertical orientation (default: FALSE)
#' @return path to an svg file
#' @importFrom fgsea plotEnrichment fgsea
#' @importFrom ggplot2 ggtitle
#' @import svglite
gseaPlot <- function(es, rankBy, selectedGenes, width, height, vertical=FALSE) {
    featureData <- fData(es)
    colnames(featureData) <- fvarLabels(es)

    ranks <- featureData[, rankBy]
    if (!is.numeric(ranks)) {
        ranks <- as.numeric(ranks)
    }
    names(ranks) <- as.character(seq_along(ranks))

    pathway <- as.character(selectedGenes)

    fgseaRes <- fgsea(list(pathway), ranks, nperm=2000, nproc=1)

    pvalString <- if (fgseaRes$nMoreExtreme == 0) {
        " < 1e-3"
    } else {
        sprintf(" = %.2g", fgseaRes$pval)
    }

    labelString <- sprintf("p-value %s, NES = %.2f", pvalString, fgseaRes$NES)


    p <- plotEnrichment(pathway, ranks) + ggtitle(labelString)
    if (vertical) {
        p <- p +
            scale_x_reverse(limits=c(length(ranks) + 1, -1), expand=c(0, 0)) +
            coord_flip() +
            NULL
    }
    p <- p + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    p
    f <- tempfile(pattern = "enrichment", tmpdir = getwd(), fileext = ".svg")
    ggsave(p, filename = f, width=width, height=height)
    # ggsave(p, filename = f, width=width, height=height)
    jsonlite::toJSON(f)
}
