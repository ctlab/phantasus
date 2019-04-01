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
#' @param quiet Boolean value which states whether the connection log should
#'     be hidden (default: TRUE)
#'
#' @return Running instance of phantasus application.
#'
#' @import opencpu
#' @import httpuv
#' @import Rook
#' @importFrom utils getFromNamespace
#' @importFrom parallel makeCluster stopCluster
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
                           openInBrowser = TRUE,
                           quiet=TRUE) {
    cacheDir <- normalizePath(cacheDir)
    preloadedDir <- if (is.null(preloadedDir))
        NULL
    else
        normalizePath(preloadedDir)

    options(phantasusCacheDir = cacheDir,
            phantasusPreloadedDir = preloadedDir)

    annotationDBMeta(cacheDir)
    FGSEAmeta(cacheDir)

    if (!opencpu:::win_or_mac()) {
        run_worker <- NULL
    } else {
        #### this fragment is adopted from opencpu::ocpu_start_server function
        #### https://github.com/opencpu/opencpu/blob/master/R/start.R
        #### :ToDo: remove code duplication

        # set root home for workers
        Sys.setenv("OCPU_MASTER_HOME" = opencpu:::tmp_root())
        on.exit(Sys.unsetenv("OCPU_MASTER_HOME"))

        # import
        sendCall <- getFromNamespace('sendCall', 'parallel')
        recvResult <- getFromNamespace('recvResult', 'parallel')
        preload <- "opencpu"

        # worker pool
        pool <- list()

        # add new workers if needed
        add_workers <- function(n = 1){
            if(length(pool) < 2){
                cl <- parallel::makeCluster(n)
                lapply(cl, sendCall, fun = function(){
                    lapply(preload, getNamespace)
                    options(phantasusCacheDir = cacheDir,
                            phantasusPreloadedDir = preloadedDir)
                    Sys.getpid()
                }, args = list())
                pool <<- c(pool, cl)
            }
        }

        # get a worker
        get_worker <- function(){
            if(!length(pool))
                add_workers(1)
            node <- pool[[1]]
            pool <<- pool[-1]
            pid <- recvResult(node)
            if(inherits(pid, "try-error"))
                warning("Worker preload error: ", pid, call. = FALSE, immediate. = TRUE)
            node$pid <- pid
            structure(list(node), class = c("SOCKcluster", "cluster"))
        }

        # main interface
        run_worker <- function(fun, ..., timeout = NULL){
            res <- tryCatch({
                if(length(timeout)){
                    setTimeLimit(elapsed = timeout)
                    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
                }
                cl <- get_worker()
                on.exit(kill_workers(cl), add = TRUE)
                node <- cl[[1]]
                sendCall(node, fun, list(...))
                recvResult(node)
            }, error = function(e){
                if(grepl("elapsed time limit", e$message)){
                    tools::pskill(node$pid)
                    stop(sprintf("Timeout reached: %ds (see rlimit.post in user.conf)", timeout))
                }
                stop(e)
            })
            if(inherits(res, "try-error"))
                stop(res)
            res
        }

        kill_workers <- function(cl){
            parallel::stopCluster(cl) # does not work when child is busy
        }

        add_workers(2)
        on.exit(kill_workers(structure(pool, class = c("SOCKcluster", "cluster"))), add = TRUE)
    }


    utils::capture.output(type = "output", {
        subPathStatic <- function (targetDirectory, subPath) {
            static <- Rook::File$new(targetDirectory)

            function (env) {
                env$PATH_INFO <- unlist(strsplit(env$PATH_INFO, subPath, fixed=TRUE))[2]
                static$call(env)
            }
        }

        app <- Rook::URLMap$new(`/phantasus/ocpu` = opencpu:::rookhandler("/phantasus/ocpu", worker_cb=run_worker),
                                `/phantasus/geo` = subPathStatic(cacheDir, '/phantasus/'),
                                `/phantasus/preloaded` = subPathStatic(cacheDir, '/phantasus/'),
                                `/phantasus/?` = subPathStatic(staticRoot, '/phantasus/'),
                                `/?` = Rook::Redirect("/phantasus/index.html"))

        tryCatch({
            server <- startServer(host, port, app = app)
            message(sprintf(
                "Server was started with following parameters: host=%s, port=%s",
                host,
                port))
        },
        error = function(e) {
            stop(paste(e,
                       "The reason may be that requested port", port,
                       "is occupied with some other application"))
        })

        if (openInBrowser) {
            url <- sprintf("http://%s:%s", host, port)
            utils::browseURL(url)
            message(paste(url, "have been opened in your default browser.\n",
                          "If nothing happened, check your 'browser'",
                          "option with getOption('browser')",
                          "or open the address manually."))
        }
        on.exit(stopServer(server))

        while(TRUE) {
            service()
            Sys.sleep(0.001)
        }
    }, split=!quiet)

}
