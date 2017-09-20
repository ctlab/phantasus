context("Load preloaded datasets")
library(Biobase)

test_that("loadPreloaded throws errors correctly", {
    expect_error(loadPreloaded("aa"), regexp = "Specify the directory")

    options(phantasusPreloadedDir = "falseDirectory")

    expect_error(loadPreloaded("aa"), regexp = "No such directory")

    options(phantasusPreloadedDir = "testdata")

    expect_error(loadPreloaded("falseFile"), regexp = "No such file")

    expect_error(loadPreloaded("wrongFormatFile"), regexp = "Wrong format")

    expect_error(loadPreloaded("wrapped_wrongFormat"), regexp = "Wrong format")

    options(phantasusPreloadedDir = NULL)
})

test_that("loadPreloaded loads file with just one ExpressionSet", {
    options(phantasusPreloadedDir = "testdata")

    expect_is(loadPreloaded("aa"), "json")

    options(phantasusPreloadedDir = NULL)
})

test_that("loadPreloaded loads file with list of ExpressionSets", {
    options(phantasusPreloadedDir = "testdata")

    expect_is(loadPreloaded("wrapped_aa"), "json")

    options(phantasusPreloadedDir = NULL)
})
