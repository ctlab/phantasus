context('counts files')
library(rhdf5)
Sys.setenv(R_USER_CONFIG_DIR = system.file("/testdata/config", package = "phantasus"))
test_that("getCountsMetaPart provides correct part of meta.rda", {
    counts_dir <- system.file("testdata/counts", package = "phantasus")
    if (counts_dir == "") {
        stop("test directory for counts does not exist")
    }
    arch_dir <- file.path(counts_dir, "archs4")
    metaPart <- phantasus:::getCountsMetaPart(counts_dir = counts_dir, collection_name = "archs4",verbose = FALSE)
    meta_file <- file.path(arch_dir, "meta.txt")
    h5_meta <- fread(meta_file, index = "file_name")
    files_in_txt <- file.path(arch_dir, h5_meta$file_name)
    files_in_meta <- file.path( counts_dir,unique(metaPart$file))
    expect_setequal(files_in_meta, files_in_txt)
    samples_in_meta <- sapply(files_in_txt, function(file) {
        sample_id <- h5_meta[file_name == basename(file)]$sample_id
        gsm_in_file <- h5read(file = file, name = sample_id)
        all(gsm_in_file %in% metaPart$accession)
    })
    expect_equal(sum(samples_in_meta), length(h5_meta$file_name))
})


test_that("updateCountsMeta works with empty directories", {
    destdir = file.path(tempdir(), "empty_dir")
    dir.create(destdir)
    updateCountsMeta(counts_dir = destdir)
    meta_file <- file.path(destdir, "meta.rda")
    expect_false(file.exists(meta_file))
    unlink(destdir,recursive = TRUE)
})
test_that("updateCountsMeta works with priorities", {
    real_path <- system.file("testdata/counts/archs4", package = "phantasus")
    if (real_path == "") {
        stop("test directory for counts does not exist")
    }
    destdir <- file.path(tempdir(), "test_counts")
    dir.create(destdir)
    tmp_path <- file.path(destdir, "archs4")
    dir.create(tmp_path)
    files <- list.files(real_path, pattern = "(.+\\.h5)|(meta.txt)", full.names = TRUE)
    file.copy(files, tmp_path)

     # without priority file
    updateCountsMeta(counts_dir = destdir)
    priority_file <- file.path(destdir, "counts_priority.txt")
    meta_file <- file.path(destdir, "meta.rda")
    expect_true(file.exists(priority_file))
    expect_true(file.exists(meta_file))
    priority <- fread(priority_file)
    expect_equal(length(priority),2)
    in_file <- priority$directory
    real_dirs <- c(".", list.dirs(destdir, full.names = FALSE, recursive =  TRUE)[-1])
    expect_true(all(in_file %in% real_dirs))

    # with prioroty file
    before_time <- file.info(priority_file)$mtime
    file.remove(meta_file)
    updateCountsMeta(counts_dir = destdir)
    expect_true(file.exists(priority_file))
    after_time <- file.info(priority_file)$mtime
    expect_equal(before_time, after_time)
    expect_true(file.exists(meta_file))

    unlink(destdir, recursive = TRUE)
})

test_that("loadCounts returns result", {
    destdir = system.file("testdata/", package = "phantasus")
    counts_dir = system.file(file.path("testdata", "counts"), package = "phantasus")
    if (counts_dir == "") {
        stop("test directory for counts does not exist")
    }
    tmp_counts <- getPhantasusConf("cache_folders")$rnaseq_counts
    if (!dir.exists(tmp_counts)){
        dir.create(tmp_counts)
    }
    dir.create(file.path(tmp_counts, "archs4"))
    count_files <- list.files(counts_dir, recursive = TRUE, full.names = FALSE)
    tmp_count_files <- file.path(tmp_counts, count_files )
    file.copy(file.path(counts_dir, count_files), tmp_count_files)
    real_wd <- getwd()
    setwd(system.file(".", package = "phantasus"))
    test_gse <- "GSE107746"
    ess <- getGSE(test_gse, mirrorPath = "https://ftp.ncbi.nlm.nih.gov")
    expect_gt(nrow(ess[[1]]), 0)
    expect_gt(ncol(ess[[1]]), 0)
    setwd(real_wd)
    unlink(tmp_counts, recursive = TRUE)

})

test_that("updateCountsMeta create good meta files for archs4", {
    destdir = system.file("testdata/", package = "phantasus")
    counts_dir = system.file(file.path("testdata", "counts"), package = "phantasus")
    if (counts_dir == "") {
        stop("test directory for counts does not exist")
    }
    tmp_dir <- file.path(tempdir(), "test_counts")
    dir.create(tmp_dir)
    tmp_counts <- file.path(tmp_dir, "counts")
    dir.create(tmp_counts)
    dir.create(file.path(tmp_counts, "archs4"))
    count_files <- list.files(counts_dir, recursive = TRUE,  pattern = ".h5$")
    tmp_count_files <- file.path(tmp_counts, count_files )
    file.copy(file.path(counts_dir, count_files), tmp_count_files)
    updateCountsMeta(counts_dir = tmp_counts)
    meta_file <- file.path(tmp_counts, "meta.rda")
    priority_file <- file.path(tmp_counts, "counts_priority.txt")
    arch_meta_file <- file.path(tmp_counts, "archs4", "meta.txt")
    expect_true(file.exists(meta_file))
    expect_true(file.exists(priority_file))
    expect_true(file.exists(arch_meta_file))
    arch_meta <- fread(arch_meta_file)
    expect_equal(nrow(arch_meta), length(count_files))
    expect_true(validateCountsCollection(file.path(tmp_counts, "archs4")))
    unlink(tmp_dir, recursive = TRUE)
})
