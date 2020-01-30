#' @param m matrix n x m
#' @param groups vector of size n of numbers from 1 to k
#' @return matrix k*m of column averages  by groups
colMeansByGroups <- function(m, groups) {
    x <- Matrix::sparseMatrix(j=seq_along(groups),
                 i=groups,
                 x=rep(1, length(groups)))
    z <- x %*% m
    res <- sweep(z, 1, table(groups), "/")
}




rasterizeHeatmap <- function(m, palette=palette, maxDimensions=c(2500, 1000)) {
    mr <- round((m - min(m))/(max(m)-min(m)) * (length(palette) - 1) + 1)
    cap <- matrix(palette[as.matrix(mr)], nrow=nrow(m))

    repN <- floor(maxDimensions/dim(m))


    cap1 <- matrix(rep(t(cap), each=repN[2]), nrow=nrow(cap), byrow = TRUE)
    cap2 <- matrix(rep(cap1, each=repN[1]), ncol=ncol(cap1), byrow = FALSE)

    # cap2[repN[1], repN[2]] == cap[1,1]
    # cap2[repN[1]+1, repN[2]] == cap[2,1]
    # cap2[repN[1], repN[2]+1] == cap[1,2]
    # cap2[repN[1]+1, repN[2]+1] == cap[2,2]

    res <- rasterGrob(cap2,
                      width=unit(1, "npc"), height=unit(1, "npc"),
                      interpolate = FALSE)
    res
}

#' Returns path to an svg file with enrichment plot
#' @param es ExpressionSet object.
#' @param rankBy name of the numeric column used for gene ranking
#' @param selectedGenes indexes of selected genes (starting from one, in the order of fData)
#' @param width width of the image (in inches)
#' @param height height of the image (in inches)
#' @param vertical whether to use vertical orientation (default: FALSE)
#' @param addHeatmap whether to add an expression heatmap, sorted by rankBy (default: FALSE)
#' @param showAnnotation a name of column annotation to add to the heatmap, default: NULL (no annotation)
#' @param pallete a vector of colors to draw heatmap
#' @param annotationColors a list of colors to use in annotation
#' @return path to an svg file
#' @importFrom fgsea plotEnrichment fgseaMultilevel
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
                     showAnnotation=NULL,
                     annotationColors=NULL,
                     pallete=c("blue", "white", "red")) {
    fullPalette <- colorRampPalette(pallete)(50)
    featureData <- fData(es)
    colnames(featureData) <- fvarLabels(es)
    absEps <- 1e-10

    ranks <- featureData[, rankBy]
    if (!is.numeric(ranks)) {
        ranks <- as.numeric(ranks)
    }
    names(ranks) <- as.character(seq_along(ranks))

    pathway <- as.character(selectedGenes)

    system.time(

    fgseaRes <- fgseaMultilevel(list(p=pathway), ranks, sampleSize = 101, nproc=1, absEps = absEps/2)
    )

    pvalString <- if (fgseaRes$pval < absEps) {
        sprintf("<%.2g", absEps)
    } else {
        sprintf("\u2248%.2g", fgseaRes$pval)
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
        aggr <- colMeansByGroups(mat, grouping)

        annotation_col <- NULL
        annotation_colors <- NULL
        if (!is.null(showAnnotation)) {
            values <- es[[showAnnotation]]
            annotation_col <- data.frame(row.names=colnames(es),
                                    setNames(list(factor(values,
                                                            levels=unique(values))),
                                             showAnnotation))

            if (!is.null(annotationColors)) {
                annotation_colors <- list()
                annotation_colors[[showAnnotation]] <- as.character(annotationColors)
                names(annotation_colors[[showAnnotation]]) <- names(annotationColors)
            }
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
                                     color=fullPalette,
                                     annotation_col = annotation_col,
                                     annotation_colors = annotation_colors,
                                     legend = FALSE,
                                     silent = TRUE)

            heatmapGrobs <- ph$gtable
            hgMatrix <- rasterizeHeatmap(aggr,
                                         palette=fullPalette,
                                         maxDimensions=c(2500, 1000))

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
                                     color=fullPalette,
                                     annotation_row = annotation_col,
                                     annotation_colors = annotation_colors,
                                     legend = FALSE,
                                     silent = TRUE)

            heatmapGrobs <- ph$gtable
            hgMatrix <- rasterizeHeatmap(Matrix::t(aggr),
                                         palette=fullPalette,
                                         maxDimensions=c(1000, 2500))

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
