context("Load preloaded datasets")
library(Biobase)

test_that("loadPreloaded throws errors correctly", {
    expect_error(loadPreloaded("aa"), regexp = "Specify the directory")

    options(phantasusPreloadedDir = "falseDirectory")
    expect_error(loadPreloaded("aa"), regexp = "No such directory")

    options(phantasusPreloadedDir = system.file("testdata", package="phantasus"))
    expect_error(loadPreloaded("falseFile"), regexp = "No such dataset")
    expect_error(loadPreloaded("wrongFormatFile"), regexp = "Wrong format")

    options(phantasusPreloadedDir = NULL)
})

test_that("loadPreloaded loads file with just one ExpressionSet", {
    options(phantasusPreloadedDir = system.file("testdata", package="phantasus"))
    expect_is(loadPreloaded("aa"), "json")

    options(phantasusPreloadedDir = NULL)
})
