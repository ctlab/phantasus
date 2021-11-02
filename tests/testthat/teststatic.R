context("Static files")

test_that("versions match", {
    versionStr <- readLines(system.file("www/phantasus.js/RELEASE.js", package="phantasus"))[1]
    versionStr <- gsub("^.*'(.*)';$", "\\1", versionStr)
    expect_equal(versionStr, fromJSON(phantasusVersion()))
})
