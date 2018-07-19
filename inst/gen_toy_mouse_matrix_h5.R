library(GEOquery)
library(rhdf5)

es <- getGEO("GSE99709")[[1]]

srcfile <- "./cache/mouse_matrix.h5"

samples <- h5read(srcfile, "meta/Sample_geo_accession")
genes <- as.character(h5read(srcfile, "meta/genes"))

keepGSMs <- head(intersect(es$geo_accession, samples), 3)

es <- es[, es$geo_accession %in% keepGSMs]

sampleIndexes <- match(es$geo_accession,
                       samples)

expression <- h5read(srcfile,
                     "data/expression",
                     index=list(seq_along(genes),
                                stats::na.omit(sampleIndexes)))
rownames(expression) <- genes
colnames(expression) <- colnames(es)[!is.na(sampleIndexes)]
H5close()


destfile <- "./inst/testdata/mouse_matrix.h5"

h5createFile(destfile)

h5createGroup(destfile, "data")

# h5createDataset(file=destfile,
#                 dataset="data/expression",
#                 dims=dim(expression),
#                 storage.mode="integer")

h5write(expression,
        file=destfile,
        name="data/expression",
        index=list(
            NULL,
            NULL))

h5createGroup(destfile, "meta")

# h5createDataset(file=destfile,
#                 dataset="meta/genes",
#                 dims=length(genes),
#                 storage.mode="character",
#                 size = max(nchar(genes)*2))

h5write(genes,
        file=destfile,
        name="meta/genes")

h5write(colnames(expression),
        file=destfile,
        name="meta/Sample_geo_accession")

H5close()
