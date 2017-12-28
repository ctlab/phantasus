#' Constructs data frame with gene annotations and submits it into Shiny GAM web-server
#' @param fData list of annotation columns
#' @param fvarLabels vector of column names
#' @return URL for Shiny GAM
#' @importFrom utils write.table
shinyGAMAnalysis <- function(fData, fvarLabels) {
    featureData <- data.frame(fData)
    colnames(featureData) <- fvarLabels
    de <- featureData


    deFile <- tempfile()
    write.table(de, deFile, sep="\t", row.names=FALSE, col.names=TRUE)

    # :ToDo: deal with bad certificates and remove ssl_verifypeer option
    r <- httr::POST(url="https://artyomovlab.wustl.edu/upload",
               body=readBin(deFile, what="raw", n=20e6),
               httr::config(ssl_verifypeer = FALSE))
    httr::stop_for_status(r)
    key <- httr::content(r, as="text", encoding="UTF-8")

    shinyGAMUrl <- sprintf("https://artyomovlab.wustl.edu/shiny/gam/?geneDE_key=%s", key)
    return(jsonlite::toJSON(shinyGAMUrl))
}
