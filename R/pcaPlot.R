#' Makes nice PCA plot for expression data
#' @param es an ExpressionSet object, should be normalized
#' @param c1 a number of the first component to plot (numeric)
#' @param c2 a number of the second component to plot (numeric)
#' @examples
#' pcaPlot(es.norm, 1, 2) + aes(color=time)
#' @export
pcaPlot <- function(es, c1, c2) {
  stopifnot(require(ggplot2))
  pca <- prcomp(t(exprs(es)))


  explained <- (pca$sdev)^2 / sum(pca$sdev^2)

  xs <- sprintf("PC%s", seq_along(explained))
  xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)

  pp <- ggplot(data=cbind(as.data.frame(pca$x), pData(es)))

  pp +
    geom_point(aes(x=eval(parse(text=xs[c1])),
                   y=eval(parse(text=xs[c2]))), size=3) +
    xlab(xlabs[c1]) + ylab(xlabs[c2])
}
