#' Serve phantasus.
#'
#' \code{servePhantasus} starts http server handling phantasus static files
#'     and opencpu server.
#'
#' @param host Host to listen.
#'
#' @param port Port to listen.
#'
#' @param staticRoot Path to static files with phantasus.js
#'     (on local file system).
#'
#' @param cacheDir Full path to cache directory.
#'
#' @param preloadedDir Full path to directory with preloaded files.
#'
#' @param runInBackground Boolean variable that states if
#'   function must run in background.
#'
#' @param openInBrowser Boolean variable that states if
#'   application should be automatically opened in default browser.
#'
#' @return Running instance of phantasus application.
#'
#' @import opencpu
#' @import httpuv
#' @import Rook
#' @export
#'
#' @examples
#' servePhantasus()
#'
servePhantasus <- function(host = '0.0.0.0', port = 8000,
                            staticRoot = system.file("www/phantasus.js",
                                                    package = "phantasus"),
                                                    cacheDir = tempdir(),
                                                    preloadedDir = NULL,
                                                    runInBackground = TRUE,
                                                    openInBrowser = TRUE) {
    options(phantasusCacheDir = cacheDir, phantasusPreloadedDir = preloadedDir)

    hostDeclaration <- paste0('host <- ', '"', host, '"', ';')
    portDeclaration <- paste0('port <- ', port, ';')
    staticRootDeclaration <- paste0('staticRoot <- ', '"', staticRoot, '"', ';')

    appDeclaration <- 'app <- Rook::URLMap$new(`/ocpu` = opencpu:::rookhandler("/ocpu"),
                            `/?` = Rook::Static$new(urls = c("/"),
                            root = staticRoot));'

    runServerScript <- 'httpuv::runServer(host, port, app = app);'

    fullScript <- c(hostDeclaration, portDeclaration,
                        staticRootDeclaration, appDeclaration, runServerScript)

    scriptFile <- tempfile(pattern = "script", fileext = ".R")

    writeLines(fullScript, scriptFile)

    if (openInBrowser) {
      browseURL(paste0(host, ':', port))
    }

    system(paste('Rscript', scriptFile), wait = !runInBackground)
}
