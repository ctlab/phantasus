context("Load preloaded datasets")
library(Biobase)

test_that("loadPreloaded throws errors correctly", {
    expect_error(loadPreloaded("aa"), regexp = "Specify the directory")

    options(phantasusPreloadedDir = "falseDirectory")

    expect_error(loadPreloaded("aa"), regexp = "No such directory")

    options(phantasusPreloadedDir = system.file("testdata", package="phantasus"))

    expect_error(loadPreloaded("falseFile"), regexp = "No such file")

    expect_error(loadPreloaded("wrongFormatFile"), regexp = "Wrong format")

    expect_error(loadPreloaded("wrapped_wrongFormat"), regexp = "Wrong format")

    options(phantasusPreloadedDir = NULL)
})

test_that("loadPreloaded loads file with just one ExpressionSet", {
    options(phantasusPreloadedDir = system.file("testdata", package="phantasus"))

    expect_is(loadPreloaded("aa"), "json")

    options(phantasusPreloadedDir = NULL)
})

test_that("loadPreloaded loads file with list of ExpressionSets", {
    options(phantasusPreloadedDir = system.file("testdata", package="phantasus"))

    expect_is(loadPreloaded("wrapped_aa"), "json")

    options(phantasusPreloadedDir = NULL)
})

test_that("checkPreloadedNames", {
    options(phantasusPreloadedDir = system.file("testdata", package="phantasus"))

    expect_equal(fromJSON(checkPreloadedNames("wrapped_aa")), c("wrapped_aa_1", "wrapped_aa_2"))
    expect_equal(fromJSON(checkPreloadedNames("aa")), c("aa"))

    expect_error(checkPreloadedNames("falseFile"), regexp = "No such file")

    options(phantasusPreloadedDir = NULL)

    expect_error(checkPreloadedNames("noDirectory"), regexp = "No such directory")

    options(phantasusPreloadedDir = "falseDirectory")

    expect_error(checkPreloadedNames("noDirectory"), regexp = "No such directory")

    options(phantasusPreloadedDir = NULL)
})
