#' @export
shinyGAMAnalysis <- function(fData, fvarLabels, orgCode) {
    featureData <- data.frame(fData)
    colnames(featureData) <- fvarLabels
    de <- featureData


    #deFile <- tempfile()
    deFile <- "/tmp/asd.txt"
    write.table(de, deFile, sep="\t", row.names=T, col.names=NA)

    # :TODO: get ssl.verifypeer from option!!!
    r <- httr::POST(url="https://artyomovlab.wustl.edu/upload",
               body=readBin(deFile, what="raw", n=20e6),
               httr::config(ssl_verifypeer = FALSE))
    httr::stop_for_status(r)
    key <- httr::content(r, as="text", encoding="UTF-8")

    shinyGAMUrl <- sprintf("https://artyomovlab.wustl.edu/shiny/gam/?oranism=%s&geneDE_key=%s", orgCode, key)
    return(jsonlite::toJSON(shinyGAMUrl))
}
