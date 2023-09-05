

#' Setup phantasus.
#' Read user config file ( or create default one) and fill \code{cache_root} using sources in \code{file}
#' @param setup_name  name of config from \code{file}. If unset or not existed, "default".
#' @param file Location of the config file
setupPhantasus <- function(
        setup_name = "default",
        file = system.file("configs/setup.yml", package = "phantasus")) {
    setup_config <-  config::get(   config = setup_name,
                                    file = file,
                                    use_parent = FALSE
                                )
    user_conf_dir <-  tools::R_user_dir(package = "phantasus", which = "config")
    user_conf_file <- file.path(user_conf_dir, "user.conf")
    cache_root = NULL
    if (!file.exists(user_conf_file)){
        cache_root = readline("Enter cache root folder for Phantasus  or leave blank for default (R tempdir):")
        if (nchar(cache_root)==0){
            cache_root = tempdir()
            }
        user_conf <- create_user_conf(cache_root = cache_root, setup_config$geo_mirrors)

        dir.create(user_conf_dir,showWarnings = FALSE, recursive = FALSE)
        cat(user_conf, file = user_conf_file)
    }
    user_conf <- getPhantasusConf()
    Sys.setenv(R_CONFIG_ACTIVE = "default")

    user_conf <- getPhantasusConf()
    message(paste("Current configuration: ", Sys.getenv("R_CONFIG_ACTIVE")))
    user_conf
    default_timeout <- getOption("timeout")
    options(timeout = max(1000, default_timeout))

    configureGEOcache(user_conf = user_conf )
    configureAnnotDB(user_conf = user_conf, setup_config = setup_config )
    configureFGSEA(user_conf = user_conf, setup_config = setup_config )
    configureRnaseqCounts(user_conf = user_conf, setup_config = setup_config )

    options(timeout =  default_timeout)
}


#' Read Phantasus Config
#'
#' @param value Value to retrieve from the config file.
#' @param configName R_CONFIG_ACTIVE value. If unset, "default".
#' @param file Location of the config file
getPhantasusConf <-  function(
        value = NULL,
        configName = Sys.getenv("R_CONFIG_ACTIVE"),
        # Modify this if your config file is somewhere else
        file = file.path(tools::R_user_dir(package = "phantasus", which = "config"), "user.conf")
) {
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
    if (dir.exists(geoPath)){
        message("! GEO cache folder exists and is treated as already configured !")
        return()
    }
    dir.create(geoPath, showWarnings = FALSE, recursive = TRUE)
}


configureAnnotDB <- function(user_conf, setup_config){
    message("Configure AnnotationDb...")
    local_path <- user_conf$cache_folders$annot_db

    if (dir.exists(local_path)){
        message("! AnnotationDB folder exists and is treated as already configured !")
        return()
    }
    menu_choices = c("Keep empty")
    actions <- c(function() {message("Kept empty")})
    dbFiles <- list.files(local_path, recursive = FALSE)
    mm_pkg <- system.file(package='org.Mm.eg.db')
    hs_pkg <- system.file(package='org.Hs.eg.db')

    if (nchar(mm_pkg) | nchar(hs_pkg)){
        message("Found local annotation databases:")
        message(paste ("\t" , mm_pkg))
        message(paste ("\t" , hs_pkg))
        menu_choices = c(menu_choices, "Use local databases")
        actions = c(actions, function(){
            if (nchar(mm_pkg)){
                library(org.Mm.eg.db)
                mm_res <- file.copy(org.Mm.eg_dbfile(), file.path(local_path, "org.Mm.eg.db"))
                if (mm_res) {
                    message(paste("! Use copy of :", org.Mm.eg_dbfile()))
                }
            }
            if (nchar(hs_pkg)){
                library(org.Hs.eg.db)
                hs_res <- file.copy(org.Hs.eg_dbfile(), file.path(local_path, "org.Hs.eg.db"))
                if (hs_res) {
                    message(paste("! Use copy of :",org.Hs.eg_dbfile()))
                }
            }
        })
    }
    if (!length(setup_config$annot_db) == 0){
        file_index <- getFileIndexDF(base_url = setup_config$annot_db, pattern = "txt|sqlite" )
        total_size <- round(sum(as.numeric(fs::fs_bytes(file_index$Size))) / 2^20, digits = 0)
        menu_choices <- c(menu_choices,  paste("Try to download from external source (~", total_size, "MB)"))
        actions <- c(actions, function(){
            loadAllFiles(url = setup_config$annot_db,
                         file_df = file_index,
                         destdir = local_path)
        } )
    }

    menu_res <- menu(choices = menu_choices, graphics = FALSE, title = "Choose how to set up the AnnotationDB folder: ")
    if (menu_res == 0){
        message("Canceled")
        return()
    }
    dir.create(local_path, recursive = TRUE)
    actions[[menu_res]]()
}


configureFGSEA <- function(user_conf, setup_config){
    message("Configure FGSEA pathways...")
    local_path <- user_conf$cache_folders$fgsea_pathways
    if (dir.exists(local_path)){
        message("! FGSEA folder exists and is treated as already configured !")
        return()
    }

    menu_choices = c("Keep empty")
    actions <- c(function() {message("Kept empty")})
    if (length(setup_config$fgsea_pathways) != 0){
        message("External source for FGSEA pathways is scpecified.")
        file_index <- getFileIndexDF(base_url = setup_config$fgsea_pathways, pattern = "txt|rds" )
        total_size <- round(sum(as.numeric(fs::fs_bytes(file_index$size))) / 2^20, digits = 2)
        menu_choices = c(menu_choices, paste0("Download files from external source? (~", total_size, "MB)"))

        actions <- c(actions, function(){
            loadAllFiles(url = setup_config$fgsea_pathways,
                         destdir = local_path,
                         file_df = file_index)
        } )
    }
    menu_res <- menu(choices = menu_choices, graphics = FALSE, title = "Choose how to set up the FGSEA pathways folder: ")
    if (menu_res == 0){
        message("Canceled")
        return()
    }
    dir.create(local_path, recursive = TRUE)
    actions[[menu_res]]()

}

configureRnaseqCounts <- function(user_conf, setup_config){
    message("Configure RNA-seq counts...")
    local_path <- user_conf$cache_folders$rnaseq_counts
    if (dir.exists(local_path)){
        message("! RNA-seq counts folder exists and is treated as already configured !")
        return()
    }
    if (length(setup_config$rnaseq_counts) == 0){
        message("! Only local RNA-seq counts will be used !")
        return()
    }

    message("External source for RNA-seq counts is scpecified.")
    menu_choices = c("Keep empty")
    actions <- c(function() {message("Kept empty")})

    req <- httr::GET(setup_config$rnaseq_counts)
    if(req$headers$server == "Highly Scalable Data Service (HSDS)"){
        message(paste(setup_config$rnaseq_counts, "response as HSDS server"))
        phantasus_lite <- system.file(package='phantasus-lite')
        if (nchar(phantasus_lite) == 0){
            message("phantasus-lite package is not installed. HSDS will be ignored")
            return()
        }
        use_hsds <- askYesNo(msg = "Would you like to use it as RNA-seq counts source?", default = FALSE)
        if (use_hsds){
            message("not implemented yet")
            #options(UseHSDS = TRUE)
        }
    } else {
        #message(paste(setup_config$rnaseq_counts, "is not HSDS server and is teated as file storage."))
        file_index <- getFileIndexDF(base_url = setup_config$rnaseq_counts, pattern = "txt|h5" )
        total_size <- round(sum(as.numeric(fs::fs_bytes(file_index$size))) / 2^20, digits = 2)
        menu_choices = c(menu_choices, paste0("Download files from external source? (~", total_size, "MB)"))

        actions <- c(actions, function(){
            loadAllFiles(url = setup_config$rnaseq_counts,
                         destdir = local_path,
                         file_df = file_index)
        } )
    }
    menu_res <- menu(choices = menu_choices, graphics = FALSE, title = "Choose how to set up the RNA-seq counts folder: ")
    if (menu_res == 0){
        message("Canceled")
        return()
    }
    dir.create(local_path, recursive = TRUE)
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
    for (file in file_df$file) {
        target_file <- file.path(destdir,file)
        dir.create(dirname(target_file), recursive = TRUE, showWarnings = FALSE)
        safeDownload(
            url = paste(url, file, sep = "/"),
            filename = basename(target_file),
            dir = dirname(target_file)
        )
    }
}


create_user_conf <- function( cache_root, geo_mirrors){

    cache_folders = list (
        geo_path = "file.path(cache_root, 'geo/')",
        rnaseq_counts = "file.path(cache_root, 'counts/')",
        annot_db  = "file.path(cache_root, 'annotationdb/')",
        fgsea_pathways = "file.path(cache_root, 'fgsea/')"
    )
    for (i in seq_along(cache_folders)){
        attr(cache_folders[[i]], "tag") <- "!expr"
    }
    static_root <- "system.file('www/phantasus.js', package = 'phantasus')"
    attr(static_root, "tag") <- "!expr"
    geo_mirrors <- if (length(geo_mirrors) == 0){
                       list(true_geo = "https://ftp.ncbi.nlm.nih.gov")
                    } else {
                       geo_mirrors
                    }
    for (i in seq_along(geo_mirrors)) {
        attr(geo_mirrors[[i]], "quoted") <- TRUE
    }
    attr(geo_mirrors, "quoted") <- TRUE
    user_conf <- list(default =
                          list(
                                host = "0.0.0.0",
                                port = "8000",
                                open_in_browser = TRUE,
                                quiet = TRUE,
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



