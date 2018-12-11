#' Returns path to an svg file with enrichment plot
#' @param es ExpressionSet object.
#' @param rankBy name of the numeric column used for gene ranking
#' @param selectedGenes indexes of selected genes (starting from one, in the order of fData)
#' @param width width of the image (in inches)
#' @param height height of the image (in inches)
#' @param vertical whether to use vertical orientation (default: FALSE)
#' @param addHeatmap whether to add an expression heatmap, sorted by rankBy (default: FALSE)
#' @param showAnnotation a name of column annotation to add to the heatmap, default: NULL (no annotation)
#' @return path to an svg file
#' @importFrom fgsea plotEnrichment fgsea
#' @importFrom ggplot2 ggtitle
#' @importFrom grDevices colorRampPalette dev.off svg
#' @importFrom utils head
#' @import grid
#' @import gtable
#' @import svglite
#' @import stats
gseaPlot <- function(es, rankBy, selectedGenes, width, height,
                     vertical=FALSE,
                     addHeatmap=FALSE,
                     showAnnotation=NULL) {
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
        "<1e-3"
    } else {
        sprintf("=%.2g", fgseaRes$pval)
    }

    labelString <- sprintf("p-value%s, NES=%.2f", pvalString, fgseaRes$NES)


    p <- plotEnrichment(pathway, ranks) + ggtitle(NULL, subtitle=labelString)
    if (vertical) {
        p <- p +
            scale_x_reverse(limits=c(length(ranks) + 1, -1), expand=c(0, 0)) +
            coord_flip() +
            NULL
    } else {
        p <- p + scale_x_continuous(limits=c(-1, length(ranks) + 1), expand=c(0, 0))
    }
    p <- p + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    enrichmentGrob <- ggplotGrob(p)

    if (addHeatmap && ncol(es) > 1) {
        # preparing heatmap
        mat <- exprs(es)[order(ranks, decreasing = TRUE), ]
        mat <- t(apply(mat, 1, scales::rescale))
        grouping <- ceiling(seq_len(nrow(mat)) / nrow(mat) * 1000)
        aggr <- Matrix.utils::aggregate.Matrix(mat, groupings=grouping, fun="mean")

        annotation_col <- NULL
        if (!is.null(showAnnotation)) {
            values <- es[[showAnnotation]]
            annotation_col <- data.frame(row.names=colnames(es),
                                    setNames(list(factor(values,
                                                            levels=unique(values))),
                                             showAnnotation))
        }


        # adding column to the left

        # 4 -> 1.5
        # 6 -> 2
        # 20+ -> 3

        heatmapWidth <- max(1, 3 - 6/ncol(aggr))

        if (vertical) {
            heatmapWidth <- min(heatmapWidth, width/4)
            ph <- pheatmap::pheatmap(aggr,
                                     cluster_rows = FALSE, cluster_cols = FALSE,
                                     show_rownames = FALSE, show_colnames = FALSE,
                                     color=colorRampPalette(c("blue", "white", "red"))(50),
                                     annotation_col = annotation_col,
                                     legend = FALSE,
                                     silent = TRUE)

            heatmapGrobs <- ph$gtable
            hgMatrix <- heatmapGrobs$grobs[[head(which(heatmapGrobs$layout$name == "matrix"), 1)]]

            panel_id <- enrichmentGrob$layout[enrichmentGrob$layout$name == "panel",c("t","l")]

            enrichmentGrob <- gtable_add_cols(enrichmentGrob, unit(heatmapWidth, "inch"), 0)
            enrichmentGrob <- gtable_add_grob(enrichmentGrob, hgMatrix,
                                              t = panel_id$t, l = 1, name="matrix")

            if (!is.null(showAnnotation)) {
                hgColAnnotation <- heatmapGrobs$grobs[[head(which(heatmapGrobs$layout$name == "col_annotation"), 1)]]
                hgAnnotationLegend <- heatmapGrobs$grobs[[head(which(heatmapGrobs$layout$name == "annotation_legend"), 1)]]
                hgLegend <- gtable_filter(heatmapGrobs, "annotation_legend", fixed=TRUE, trim=TRUE)

                enrichmentGrob <- gtable_add_grob(enrichmentGrob, hgColAnnotation,
                                                  t = 1, l = 1, b=panel_id$t-1, name="col_annotation")


                enrichmentGrob <- gtable_add_cols(enrichmentGrob, gtable_width(hgLegend), 0)
                enrichmentGrob <- gtable_add_grob(enrichmentGrob, hgAnnotationLegend,
                                                  t = panel_id$t, l = 1, name="annotation_legend")

            }

        } else {
            # horizontal
            heatmapWidth <- min(heatmapWidth, height/4)
            ph <- pheatmap::pheatmap(Matrix::t(aggr),
                                     cluster_rows = FALSE, cluster_cols = FALSE,
                                     show_rownames = FALSE, show_colnames = FALSE,
                                     color=colorRampPalette(c("blue", "white", "red"))(50),
                                     annotation_row = annotation_col,
                                     legend = FALSE,
                                     silent = TRUE)

            heatmapGrobs <- ph$gtable
            hgMatrix <- heatmapGrobs$grobs[[head(which(heatmapGrobs$layout$name == "matrix"), 1)]]

            panel_id <- enrichmentGrob$layout[enrichmentGrob$layout$name == "panel",c("t","l")]


            enrichmentGrob <- gtable_add_rows(enrichmentGrob, unit(heatmapWidth, "inch"), 0)
            enrichmentGrob <- gtable_add_grob(enrichmentGrob, hgMatrix,
                                              t = 1, l=panel_id$l, name="matrix")

            if (!is.null(showAnnotation)) {
                hgColAnnotation <- heatmapGrobs$grobs[[head(which(heatmapGrobs$layout$name == "row_annotation"), 1)]]
                hgAnnotationLegend <- heatmapGrobs$grobs[[head(which(heatmapGrobs$layout$name == "annotation_legend"), 1)]]
                hgLegend <- gtable_filter(heatmapGrobs, "annotation_legend", fixed=TRUE, trim=TRUE)

                enrichmentGrob <- gtable_add_grob(enrichmentGrob, hgColAnnotation,
                                                  t = 1, l = 1, b=1, r=panel_id$l-1, name="col_annotation")


                enrichmentGrob <- gtable_add_cols(enrichmentGrob, gtable_width(hgLegend), 0)
                enrichmentGrob <- gtable_add_grob(enrichmentGrob, hgAnnotationLegend,
                                                  t = 1, l = 1, name="annotation_legend")

            }

        }





    }

    f <- tempfile(pattern = "enrichment", tmpdir = getwd(), fileext = ".svg")
    svg(f, width=width, height=height)
    grid.draw(enrichmentGrob)
    dev.off()
    # ggsave(p, filename = f, width=width, height=height)
    jsonlite::toJSON(f)
}
