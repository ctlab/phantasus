#' Constructs data frame with gene annotations and submits it into Shiny GAM web-server
#' @return URL for Shiny GAM
#' @import httr
shinyGAMAnalysis <- function(fData, fvarLabels) {
    featureData <- data.frame(fData)
    colnames(featureData) <- fvarLabels
    de <- featureData


    deFile <- tempfile()
    write.table(de, deFile, sep="\t", row.names=F, col.names=T)

    # :ToDo: deal with bad certificates and remove ssl_verifypeer option
    r <- httr::POST(url="https://artyomovlab.wustl.edu/upload",
               body=readBin(deFile, what="raw", n=20e6),
               httr::config(ssl_verifypeer = FALSE))
    httr::stop_for_status(r)
    key <- httr::content(r, as="text", encoding="UTF-8")

    shinyGAMUrl <- sprintf("https://artyomovlab.wustl.edu/shiny/gam/?geneDE_key=%s", key)
    return(jsonlite::toJSON(shinyGAMUrl))
}
