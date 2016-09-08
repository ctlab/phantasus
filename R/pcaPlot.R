#' Makes nice PCA plot for expression data
#' @param es an ExpressionSet object, should be normalized
#' @param c1 a number of the first component to plot (numeric)
#' @param c2 a number of the second component to plot (numeric)
#' @examples
#' pcaPlot(es.norm, 1, 2) + aes(color=time)
#' @export
pcaPlot <- function(es, columns=c(), rows=c(), c1, c2, size="", colour="", label="") {
  n1 <- as.numeric(c1)
  n2 <- as.numeric(c2)
  stopifnot(require(ggplot2))
  stopifnot(require(ggrepel))
  stopifnot(require(Biobase))
  stopifnot(require(svglite))
  
  if (is.null(rows)) {
    rows <- 1:nrow(exprs(es))
  }
  if (is.null(columns)) {
    columns <- 1:ncol(exprs(es))
  }
  
  rows <- as.numeric(rows)
  columns <- as.numeric(columns)
  
  data <- t(exprs(es)[rows,columns])
  pca <- prcomp(data)
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)
  
  xs <- sprintf("PC%s", seq_along(explained))
  xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)
  
  pData <- pData(es)[!(rownames(pData(es)) %in% setdiff(rownames(pData(es)), rownames(pca$x))),]
  
  if (size != "") {
    pData[[size]] <- as.numeric(pData[[size]])
  }
  
  pp <- ggplot(data=cbind(as.data.frame(pca$x), pData, sampleNames(es)))
  if (size != "" && colour != "") {
    
    aes <- aes_string(x=xs[n1],
                      y=xs[n2], colour=colour, size=size)
    
  } else if (colour != "") {
    aes <- aes_string(x=xs[n1],
                      y=xs[n2], colour=colour)
  } else if (size != "") {
    aes <- aes_string(x=xs[n1],
                      y=xs[n2], size=size)
  } else {
    aes <- aes_string(x=xs[n1],
                      y=xs[n2])
  }
  
  g <- pp + aes + 
    geom_point() + 
    xlab(xlabs[n1]) + ylab(xlabs[n2])
  
  if (label == "id") {
    label <- "sampleNames(es)"
  }
  if (label != "") {    
    message("i'm here 2s")
    g <- g + geom_text_repel(aes_string(label=label))
    
  }
  f <- tempfile(pattern="plot",tmpdir=getwd(),fileext=".svg")
  ggsave(f, g)
  
  
  
  print(capture.output(str(g)))
  return(f)
}
