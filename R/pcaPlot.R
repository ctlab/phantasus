#' Makes nice PCA plot for expression data
#' @param es an ExpressionSet object, should be normalized
#' @param c1 a number of the first component to plot (numeric)
#' @param c2 a number of the second component to plot (numeric)
#' @examples
#' pcaPlot(es.norm, 1, 2) + aes(color=time)
#' @export
pcaPlot <- function(es, columns, c1, c2, size="", colour="") {
  stopifnot(require(ggplot2))
  stopifnot(require(Biobase))
  if (is.na(size)) {
    size = ""
  }
  if (is.na(colour)) {
    colour = ""
  }
  data <- t(exprs(es)[columns,])
  pca <- prcomp(~., data.frame(data))
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)
  print(size)
  print(colour)

  xs <- sprintf("PC%s", seq_along(explained))
  xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)

  pData <- pData(es)[!(rownames(pData(es)) %in% setdiff(rownames(pData(es)), rownames(pca$x))),]
  pp <- ggplot(data=cbind(as.data.frame(pca$x), pData))
  if (size != "" && colour != "") {
    
    pData[[size]] <- as.numeric(pData[[size]])
    
    aes <- aes_string(x=xs[c1],
                      y=xs[c2], colour=colour, size=size)
    
  } else if (colour != "") {
    aes <- aes_string(x=xs[c1],
               y=xs[c2], colour=colour)
  } else if (size != "") {
    pData[[size]] <- as.numeric(pData[[size]])
    aes <- aes_string(x=xs[c1],
        y=xs[c2], size=size)
  } else {
    aes <- aes_string(x=xs[c1],
                      y=xs[c2])
  }
  
  
  pp +
    geom_point(aes) +
    xlab(xlabs[c1]) + ylab(xlabs[c2])
}
