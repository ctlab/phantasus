context("GSEA plot")
library(jsonlite)
library(Biobase)

test_that("gseaPlot works", {
    load(file = system.file("testdata/GSE27112-GPL6103.rda", package="phantasus"))

    set.seed(42)
    fData(es)$t <- rnorm(nrow(es))

    f <- fromJSON(gseaPlot(es, rankBy = "t", selectedGenes = sample.int(nrow(es), 10),
             width=6, height=4))
    expect_true(file.exists(f))

    f <- fromJSON(gseaPlot(es, rankBy = "t", selectedGenes = sample.int(nrow(es), 10),
             width=6, height=4, vertical = TRUE))
    expect_true(file.exists(f))

    f <- fromJSON(gseaPlot(es, rankBy = "t", selectedGenes = sample.int(nrow(es), 10),
             width=6, height=4, vertical = TRUE,
             addHeatmap = TRUE, showAnnotation = "time"))
    expect_true(file.exists(f))

    f <- fromJSON(gseaPlot(es, rankBy = "t", selectedGenes = sample.int(nrow(es), 10),
             width=6, height=4, vertical = FALSE,
             addHeatmap = TRUE, showAnnotation = "time"))
    expect_true(file.exists(f))
})

