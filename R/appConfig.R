

#' Setup phantasus.
#' Read user config file ( or create default one) and fill \code{cache_root} using sources in \code{file}.
#' @param setup_name  name of config from \code{file}. If unset or not existed, "default".
#' @param file Location of the setup.yml file with setup parameters. If not existed use file from package
#' @export
setupPhantasus <- function(
        setup_name = "default",
        file = confFile("setup.yml")){
    setup_config <- getSetupConf(file, setup_name)
    user_conf_file <- confFile("user.conf")
    cache_root = ""
    if (!file.exists(user_conf_file)){
        default_cache <- tools::R_user_dir("phantasus", which = "cache")
        if (interactive()){
            message("Specify the cache root for Phantasus. \n This folder will contain all Phantasus-related files: \n - GEO cache files \n - Annotation databases \n - FGSEA pathways \n - RNA-seq counts")
            cache_root <- readline(paste0("Enter path or leave blank for default (", default_cache, "):"))
        }
        if (is.null(cache_root) || is.na(cache_root) || (nchar(cache_root)==0)){
            cache_root <- default_cache
        }
        user_conf <- get_user_conf(cache_root = cache_root, setup_config)
        dir.create(dirname(user_conf_file), showWarnings = FALSE, recursive = TRUE)
        message(paste("Create user configuration file:", user_conf_file))
        cat(user_conf, file = user_conf_file)

    } else{
        message(paste("Use existed configuration file:", user_conf_file))
    }
    user_conf <- getPhantasusConf()
    Sys.setenv(R_CONFIG_ACTIVE = "default")

    user_conf <- getPhantasusConf()
    message(paste("Current configuration: ", Sys.getenv("R_CONFIG_ACTIVE")))
    default_timeout <- getOption("timeout")
    options(timeout = max(1000, default_timeout))

    configureGEOcache(user_conf = user_conf )
    configureAnnotDB(user_conf = user_conf, setup_config = setup_config )
    configureFGSEA(user_conf = user_conf, setup_config = setup_config )
    configureRnaseqCounts(user_conf = user_conf, setup_config = setup_config )

    options(timeout =  default_timeout)
    if (file.exists(user_conf_file)){
        message(paste("! Successfully setup configuration file:", user_conf_file, "!"))
    }
}

confFile <- function(name){
    return(file.path(tools::R_user_dir(package = "phantasus", which = "config"), name ))

}

getSetupConf <- function(file, setup_name){
    if (!file.exists(file)){
        message("Use default setup parameters")
        file <- system.file("configs/R/phantasus/setup.yml", package = "phantasus")
    }
    setup_config <-  config::get(   config = setup_name,
                                    file = file,
                                    use_parent = FALSE
    )
    return(setup_config)
}

#' Read Phantasus Config
#'
#' @param value Value to retrieve from the config file.
#' @param configName R_CONFIG_ACTIVE value. If unset, "default".
#' @param file Location of the config file
#' @export
getPhantasusConf <-  function(
        value = NULL,
        configName = Sys.getenv("R_CONFIG_ACTIVE"),
        # Modify this if your config file is somewhere else
        file = file.path(tools::R_user_dir(package = "phantasus", which = "config"), "user.conf")
) {
    if (!file.exists(file)){
        stopPhantasus()
    }
    config::get(
        value = value,
        config = configName,
        file = file,
        use_parent = FALSE
    )

}

configureGEOcache <- function(user_conf){
    message("Configure GEO cache folder...")
    geoPath <- user_conf$cache_folders$geo_path
    if (dir.exists(geoPath) && rw_dir_check(geoPath)){
        message("! GEO cache folder exists and is treated as already configured !")
        return()
    }
    unsafe_dir_create(geoPath)
}


configureAnnotDB <- function(user_conf, setup_config){
    message("Configure AnnotationDb...")
    local_path <- user_conf$cache_folders$annot_db

    if (dir.exists(local_path) && rw_dir_check(local_path)){
        message("! AnnotationDB folder exists and is treated as already configured !")
        return()
    }
    menu_choices = c()
    actions <- c()
    dbFiles <- list.files(local_path, recursive = FALSE)
    mm_pkg <- system.file(package='org.Mm.eg.db')
    hs_pkg <- system.file(package='org.Hs.eg.db')

    if (nchar(mm_pkg) | nchar(hs_pkg)){
        message("Found local annotation databases:")
        message(paste ("\t" , mm_pkg))
        message(paste ("\t" , hs_pkg))
        menu_choices = c(menu_choices, "Use local databases listed above")
        actions = c(actions, function(){
            if (nchar(mm_pkg)){
                library(org.Mm.eg.db)
                mm_res <- file.copy(org.Mm.eg_dbfile(), file.path(local_path, "org.Mm.eg.sqlite"))
                if (mm_res) {
                    message(paste("! Use copy of :", org.Mm.eg_dbfile()))
                }
            }
            if (nchar(hs_pkg)){
                library(org.Hs.eg.db)
                hs_res <- file.copy(org.Hs.eg_dbfile(), file.path(local_path, "org.Hs.eg.sqlite"))
                if (hs_res) {
                    message(paste("! Use copy of :",org.Hs.eg_dbfile()))
                }
            }
        })
    }
    if (!length(setup_config$annot_db) == 0){
        file_index <- getFileIndexDF(base_url = setup_config$annot_db, pattern = "txt|sqlite" )
        total_size <- round(sum(as.numeric(fs::fs_bytes(file_index$size))) / 2^20, digits = 0)
        menu_choices <- c(menu_choices,  paste("Download files from the ctlab Phantasus mirror (~", total_size, "MB)"))
        actions <- c(actions, function(){
            loadAllFiles(url = setup_config$annot_db,
                         file_df = file_index,
                         destdir = local_path)
        } )
    }
    menu_choices = c(menu_choices, "Keep empty")
    actions <- c(actions, function() {message("Kept empty")})
    if ((!interactive()) | length(menu_choices) == 1){
        menu_res <- 1
    } else{
        menu_res <- utils::menu(choices = menu_choices, graphics = FALSE, title = "Choose how to set up the AnnotationDB folder: ")
        if (menu_res == 0){
            message("Canceled")
            return()
        }
    }
    unsafe_dir_create(local_path)
    actions[[menu_res]]()
    annotationDBMeta(getPhantasusConf("cache_folders")$annot_db)
}


configureFGSEA <- function(user_conf, setup_config){
    message("Configure FGSEA pathways...")
    local_path <- user_conf$cache_folders$fgsea_pathways
    if (dir.exists(local_path) && rw_dir_check(local_path)){
        message("! FGSEA folder exists and is treated as already configured !")
        return()
    }

    menu_choices = c()
    actions <- c()
    if (length(setup_config$fgsea_pathways) != 0){
        message("External source for FGSEA pathways is scpecified.")
        file_index <- getFileIndexDF(base_url = setup_config$fgsea_pathways, pattern = "txt|rds" )
        total_size <- round(sum(as.numeric(fs::fs_bytes(file_index$size))) / 2^20, digits = 2)
        menu_choices <- c(menu_choices, paste0("Download files from the ctlab Phantasus mirror (~", total_size, "MB)"))

        actions <- c(actions, function(){
            loadAllFiles(url = setup_config$fgsea_pathways,
                         destdir = local_path,
                         file_df = file_index)
        } )
    }
    menu_choices = c(menu_choices,"Keep empty")
    actions <- c(actions, function() {message("Kept empty")})
    if ((!interactive()) | length(menu_choices) == 1){
        menu_res <- 1
    } else{
        menu_res <- utils::menu(choices = menu_choices, graphics = FALSE, title = "Choose how to set up the FGSEA pathways folder: ")
        if (menu_res == 0){
            message("Canceled")
            return()
        }
    }
    unsafe_dir_create(local_path)
    actions[[menu_res]]()
    FGSEAmeta(getPhantasusConf("cache_folders")$fgsea_pathways)
}

configureRnaseqCounts <- function(user_conf, setup_config){
    message("Configure RNA-seq counts...")
    selected_path <- user_conf$cache_folders$rnaseq_counts
    options(PhantasusUseHSD = NULL)
    if (is.null(selected_path) || is.na(selected_path) ||  (nchar(selected_path) == 0)){
        stop("!Configuration failed: RNA-seq counts path should be specified !")
        return()
    }
    if (dir.exists(selected_path) && rw_dir_check(selected_path)){
        message("! RNA-seq counts folder exists. Update meta... !")
        updateCountsMeta(counts_dir =  selected_path, force = FALSE, verbose = FALSE)
        message("! RNA-seq counts configured !")
        return()
    }

    if(!grepl(pattern = "^http(s)?://", x = selected_path)){
        message("! RNA-seq count path looks like local path !")
        unsafe_dir_create(selected_path)
        message("RNA-seq counts configured")
        return()
    }
    message("! RNA-seq counts path looks like web URL !")
    if( !isHSDS(selected_path)){
        message("Resource doesn't respond as HSDS server")
        stop("Configuration failed: RNA-seq count path should be local diractory or HSDS server")
    }
    message(paste(selected_path, "response as HSDS server"))

    menu_choices = c("Yes, load expression matrix from remote server when the known RNA-seq dataset is requested")
    actions <- c(function() {
        phantasus_lite <- system.file(package='phantasusLite')
        if (nchar(phantasus_lite) == 0){
            message("phantasus-lite package is not installed. HSDS will be ignored")
            stop("Configuration failed: Install phantasusLite package to use HSDS server")
        }
        options(PhantasusUseHSDS = TRUE)
        message("HSDS server will be used as source of RNA-seq counts")
        })

    menu_choices <- c(menu_choices, "No, keep all RNA-seq datasets without expression data")
    actions <- c(actions, function() {
        options(PhantasusUseHSDS = FALSE)
        message("HSDS server will be ignored. RNA-seq datasets will be loaded without expression matrices.")
    })
    if (!interactive()){
        menu_res <- 1
    } else {
        menu_res <- utils::menu(choices = menu_choices, graphics = FALSE, title = "Current Phantasus configuration allows it to load count matrices from remote server when RNA-seq dataset is requested.\nWould you like to use this feature?")
        if (menu_res == 0){
            message("Canceled")
            return()
        }
    }
    actions[[menu_res]]()

}

getFileIndexDF <- function(base_url, pattern, location = "/"){
    req <- httr::GET(paste0(base_url, location))
    text <- httr::content(req, "text", encoding = "UTF-8")
    index_df <-  data.table(na.omit(XML::readHTMLTable(text, skip = c(1,2), as.data.frame = TRUE, trim = TRUE )[[1]]))
    file_df <- index_df[grepl(pattern = pattern, x = index_df$Name) & index_df$Size != "-", .(file = paste0(location, Name), date = `Last modified`, size = Size)]

    for (dir in index_df[index_df$Size == "-",]$Name){
        local_dir <- paste0(location, dir)
        file_df <- rbind(file_df, getFileIndexDF(base_url, pattern, local_dir ))
    }
    return(file_df)
}

loadAllFiles <- function(url, file_df, destdir ){
    message(paste("Trying download from", url))
    url <- trimws(url, which = "right", whitespace = "/")
    for (file in file_df$file) {
        target_file <- file.path(destdir,file)
        dir.create(dirname(target_file), recursive = TRUE, showWarnings = FALSE)
        safeDownload(
            url = paste0(url, file),
            filename = basename(target_file),
            dir = dirname(target_file)
        )
    }
}


get_user_conf <- function( cache_root, setup_config){

    cache_folders = list (
        geo_path = "file.path(cache_root, 'geo/')",
        annot_db  = "file.path(cache_root, 'annotationdb/')",
        fgsea_pathways = "file.path(cache_root, 'fgsea/')"
    )
    for (i in seq_along(cache_folders)){
        attr(cache_folders[[i]], "tag") <- "!expr"
    }
    if ( is.null(setup_config$rnaseq_counts) | length(setup_config$rnaseq_counts) == 0){
        cache_folders$rnaseq_counts <- "file.path(cache_root, 'counts/')"
        attr(cache_folders$rnaseq_counts, "tag") <- "!expr"
    } else{
        cache_folders$rnaseq_counts <- setup_config$rnaseq_counts
    }

    static_root <- "system.file('www/phantasus.js', package = 'phantasus')"
    attr(static_root, "tag") <- "!expr"
    geo_mirrors <- if (length(setup_config$geo_mirrors) == 0){
                       list(true_geo = "https://ftp.ncbi.nlm.nih.gov")
                    } else {
                        setup_config$geo_mirrors
                    }
    for (i in seq_along(geo_mirrors)) {
        attr(geo_mirrors[[i]], "quoted") <- TRUE
    }
    attr(geo_mirrors, "quoted") <- TRUE
    user_conf <- list(default =
                          list(
                                host = "0.0.0.0",
                                port = "8000",
                                preloaded_dir = NULL,
                                static_root = static_root,
                                cache_root = cache_root,
                                cache_folders = cache_folders,
                                geo_mirrors = geo_mirrors
                            )
                    )
    attr(user_conf$default$cache_root, "quoted") <- TRUE
    attr(user_conf$default$host,"quoted") <- TRUE
    class(user_conf$default$port) = 'verbatim'
    yaml::as.yaml(user_conf, handlers = list(
        logical = function(x){
            result <- ifelse(x, "TRUE", "FALSE")
            class(result) <- "verbatim"
            return(result)
        },
        "NULL" = function(x){
            result = "NULL"
            class(result) <- "verbatim"
            return(result)
        }
    ) )
}


#' Creates default docker conf file
#' Function creates default docker user configuration file based on provided \code{setup_file}
#' or on default parameters if \code{setup_file} doesn't exist. If \code{user_conf_file} exists function does nothing.
#' @param setup_file  name of config from \code{file}. If unset or not existed, "default".
#' @param user_conf_file Location of the setup.yml file with setup parameters. If not existed use file from package
createDockerConf <- function( setup_file = confFile("setup.yaml"), user_conf_file = confFile("user.conf")){
    if (!file.exists(user_conf_file)){
        setup_conf <-  getSetupConf(setup_file, "default")
        user_conf <- get_user_conf("/var/phantasus/cache", setup_conf)
        user_conf$preloaded_dir <- "/var/phantasus/preloaded"
        dir.create(dirname(user_conf_file), showWarnings = FALSE, recursive = TRUE)
        message(paste("Create user configuration file:", user_conf_file))
        cat(user_conf, file = user_conf_file)
    }
}

unsafe_dir_create <- function(dir_name,  recursive = TRUE){
    configured <- dir.create(dir_name, recursive = recursive, showWarnings = FALSE)
    if (!configured){
        stop(paste("Configuration failed: can't create", dir_name))
    }
    rw_dir_check(dir_name)
}

