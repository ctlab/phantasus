#' Starts http server handling morpheus static files and opencpu
#' @param host host to listen
#' @param port port to listen
#' @param staticRoot path to static files with phantasus.js (on local file system)
#' @param cacheDir full path to cache directory
#' @import opencpu
#' @import httpuv
#' @import Rook
#' @export
#' @examples
#' servePhantasus("0.0.0.0", 8000, cacheDir=file.path(getwd(), "cache"))
servePhantasus <- function(host, port,
                          staticRoot=system.file("www/phantasus.js", package="phantasus"),
                          cacheDir=tempdir()) {
  options(phantasusCacheDir=cacheDir)
  app <-
    Rook::URLMap$new(
      "/ocpu"=opencpu:::rookhandler("/ocpu"),
      "/?"=Rook::Static$new(
        urls = c('/'),
        root = staticRoot
      ))

  httpuv::runServer(host,
                    port,
                    app=app)

}