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
#' @param openInBrowser Boolean value which states if application will be
#'     automatically loaded in default browser.
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
#' servePhantasus()
#' }
servePhantasus <- function(host = '0.0.0.0',
                           port = 8000,
                           staticRoot = system.file("www/phantasus.js",
                                                    package = "phantasus"),
                           cacheDir = tempdir(),
                           preloadedDir = NULL,
                           openInBrowser = TRUE) {
    options(phantasusCacheDir = normalizePath(cacheDir),
            phantasusPreloadedDir = if (is.null(preloadedDir)) NULL else normalizePath(preloadedDir))


    utils::capture.output(type = "output", {
        app <- Rook::URLMap$new(`/ocpu` = opencpu:::rookhandler("/ocpu"),
                                `/?` = Rook::Static$new(urls = c("/"),
                                                        root = staticRoot))

        tryCatch({
            server <- startServer(host, port, app = app)
            message(sprintf("Server was started with following parameters: host=%s, port=%s", host, port))
        },
        error = function(e) {
            stop(paste(e,
                       "The reason may be that requested port", port, "is occupied with some other application"))
        })

        if (openInBrowser) {
            url <- sprintf("http://%s:%s", host, port)
            utils::browseURL(url)
            message(paste(url, "have been opened in your default browser.\n",
                          "If nothing happened, check your 'browser' option with getOption('browser')",
                          "or open the address manually."))
        }
        on.exit(stopServer(server))

        while(TRUE) {
            service()
            Sys.sleep(0.001)
        }
    })

}
