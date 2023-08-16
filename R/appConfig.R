

#' apply config file
#'
#' @param configName R_CONFIG_ACTIVE value. If unset, "default".
#' @param file Location of the config file
configurePhantasus <- function(
        configName = Sys.getenv("R_CONFIG_ACTIVE"),
        # Modify this if your config file is somewhere else
        file = system.file("configs/phantasus-conf.yml", package = "phantasus")
) {
    if (configName == ""){
        Sys.setenv(R_CONFIG_ACTIVE = "default")
    }
    if (configName !=  Sys.getenv("R_CONFIG_ACTIVE")){
        Sys.setenv(R_CONFIG_ACTIVE = configName)
    }
    curConf <- getPhantasusConf(file = file)
    if (!dir.exists(curConf$local_cache$cache_root)){
        dir.create(curConf$local_cache$cache_root)
    }
    configureAnnotDB(curConf)
    configureFGSEA(curConf)
    configureRnaseqCounts(curConf)
    message(paste("curent configuration: ", Sys.getenv("R_CONFIG_ACTIVE")))
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
        file = system.file("configs/phantasus-conf.yml", package = "phantasus")
) {
    config::get(
        value = value,
        config = configName,
        file = file,
        use_parent = FALSE
    )

}


configureAnnotDB <- function(curConf){
    message("Configure AnnotationDb...")
    local_path <- curConf$local_cache$annot_db

    if (length(local_path) == 0){
        local_path <- file.path(curConf$local_cache$cache_root, "annotationdb")
    }
    if (!dir.exists(local_path)){
        dir.create(local_path)
    }
    dbFiles <- list.files(local_path, recursive = FALSE)
    if (length(dbFiles) > 0) {
        message("! Local annotationDB folder is not empty and is treated as already configured !")
        return()
    }
    mm_pkg <- system.file(package='org.Mm.eg.db')
    hs_pkg <- system.file(package='org.Hs.eg.db')
    if (length(curConf$external_sources$annot_db) == 0){
        if (nchar(mm_pkg)){
            library(org.Mm.eg.db)
            mm_res <- file.symlink(org.Mm.eg_dbfile(), file.path(local_path, "org.Mm.eg.db"))
            if (mm_res) {
                message(paste("! Use :", org.Mm.eg_dbfile()))
            }
        }
        if (nchar(hs_pkg)){
            library(org.Hs.eg.db)
            hs_res <- file.symlink(org.Hs.eg_dbfile(), file.path(local_path, "org.Hs.eg.db"))
            if (hs_res) {
                message(paste("! Use :",org.Hs.eg_dbfile()))
            }
        }
    } else {
        message("External source for AnnotationDB is scpecified. Local packages will be ignored.")
        msg <- "Would you like to download files from external source (may take a while)?"
        tryLoad <- askYesNo(msg = msg, default = FALSE)
        if (tryLoad) {
            loadAllFiles( url = curConf$external_sources$annot_db,
                          destdir = local_path ,
                          pattern = "txt|sqlite")
        }
    }
}


configureFGSEA <- function(curConf){
    message("Configure FGSEA pathways...")
    local_path <- curConf$local_cache$fgsea_pathways
    if (length(local_path) == 0){
        local_path <- file.path(curConf$local_cache$cache_root, "fgsea")
    }
    if (!dir.exists(local_path)){
        dir.create(local_path)
    }
    fgseaFiles <- list.files(local_path, recursive = FALSE)
    if (length(fgseaFiles) > 0) {
        message("! Local fgsea folder is not empty and is treated as already configured !")
        return()
    }
    if (length(curConf$external_sources$fgsea_pathways) != 0){
        message("External source for FGSEA pathways is scpecified.")

        loadAllFiles(url = curConf$external_sources$fgsea_pathways,
                     destdir = local_path,
                     pattern = "txt|rds")
    }
}

configureRnaseqCounts <- function(curConf){
    message("Configure RNA-seq counts...")
    local_path <- curConf$local_cache$rnaseq_counts
    if (length(local_path) == 0){
        local_path <- file.path(curConf$local_cache$cache_root, "counts")
    }
    if (!dir.exists(local_path)){
        dir.create(local_path)
    }
    if (length(curConf$external_sources$rnaseq_counts) == 0){
        message("! Only local RNA-seq counts will be used !")
        return()
    }
    message("External source for RNA-seq counts is scpecified.")
    req <- httr::GET(curConf$external_sources$rnaseq_counts)
    if(req$headers$server == "Highly Scalable Data Service (HSDS)"){
        message(paste(curConf$external_sources$rnaseq_counts, "response as HSDS server"))
        phantasus_lite <- system.file(package='phantasus-lite')
        if (nchar(phantasus_lite) == 0){
            message("phantasus-lite package is not installed. HSDS will be ignored")
            return()
        }
        use_hsds <- askYesNo(msg = "Would you like to use it as RNA-seq counts source?", default = FALSE)
        if (use_hsds){
            options(UseHSDS = TRUE)
        }
    } else {
        message(paste(curConf$external_sources$rnaseq_counts, "is not HSDS server and is teated as file storage."))
        count_files <- list.files(local_path, recursive = TRUE,
                                  pattern = '\\.h5$')
        if(length(count_files)){
            message("! Local RNA-seq counts folder is not empty and is treated as already configured !")
            return()
        }
        message("! ATTENTION: The size of RNA-seq count files may exeed 100 GB in total ! ")
        tryLoad <- askYesNo(msg = "Would you like to download files from external source?", default = FALSE)
        if (tryLoad) {
            loadAllFiles( url = curConf$external_sources$rnaseq_counts,
                          destdir = local_path ,
                          pattern = "txt|h5")
        }
    }
}

loadAllFiles <- function(url, destdir, pattern = ".*"){
    message(paste("Trying download all from ", url))
    if (!dir.exists(destdir)){
        dir.create(destdir)
    }
    req <- httr::GET(url)
    text <- httr::content(req, "text", encoding = "UTF-8")
    index_df <-  na.omit(XML::readHTMLTable(text, skip = c(1,2), as.data.frame = TRUE, trim = TRUE )[[1]])
    file_df <- index_df[grepl(pattern = pattern, x = index_df$Name) & index_df$Size != "-",]

    for (file in file_df$Name) {
        safeDownload(
            url = paste(url, file, sep = "/"),
            filename = file,
            dir = destdir
        )
    }
    for (dir in index_df[index_df$Size == "-",]$Name){
        local_dir <- file.path(destdir, dir)
        loadAllFiles(url = paste(url, dir, sep ="/"),
                     destdir = local_dir,
                     pattern = pattern )
    }
}




