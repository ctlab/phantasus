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
#' @return Running instance of phantasus application.
#'
#' @import opencpu
#' @import httpuv
#' @import Rook
#' @export
#'
#' @examples
#' \dontrun{
#' servePhantasus('0.0.0.0', 8000, cacheDir=file.path(getwd(), 'cache'))
#' }
servePhantasus <- function(host, port,
                            staticRoot = system.file("www/phantasus.js",
                                                    package = "phantasus"),
    cacheDir = tempdir()) {
    options(phantasusCacheDir = cacheDir)
    app <- Rook::URLMap$new(`/ocpu` = opencpu:::rookhandler("/ocpu"),
                            `/?` = Rook::Static$new(urls = c("/"),
                            root = staticRoot))

    httpuv::runServer(host, port, app = app)

}
