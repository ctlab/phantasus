#' Starts http server handling morpheus static files and opencpu
#' @param host host to listen
#' @param port port to listen
#' @param morpheusRoot path to static files with morpheus (on local file system)
#' @import opencpu
#' @import httpuv
#' @import Rook
#' @export
serveMorpheus <- function(host, port, morpheusRoot) {
  app <-
    Rook::URLMap$new(
      "/ocpu"=opencpu:::rookhandler("/ocpu"),
      "/?"=Rook::Static$new(
        urls = c('/'),
        root = morpheusRoot
      ))

  httpuv::runServer(host,
                    port,
                    app=app)

}
