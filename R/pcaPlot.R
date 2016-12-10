#' Makes nice PCA plot for expression data
#' @param es an ExpressionSet object, should be normalized
#' @param c1 a number of the first component to plot (numeric)
#' @param c2 a number of the second component to plot (numeric)
#' @examples
#' pcaPlot(es.norm, 1, 2) + aes(color=time)
#' @export
pcaPlot <- function(es, columns=c(), rows=c(), c1, c2, size="", colour="", label="", replacena = "mean") {
  n1 <- as.numeric(c1)
  n2 <- as.numeric(c2)
  stopifnot(require(ggplot2))
  stopifnot(require(ggrepel))
  stopifnot(require(Biobase))
  stopifnot(require(svglite))
  stopifnot(require(plotly))
  stopifnot(require(htmltools))
  if (is.null(rows)) {
    rows <- 1:nrow(exprs(es))
  }
  if (is.null(columns)) {
    columns <- 1:ncol(exprs(es))
  }

  rows <- as.numeric(rows)
  columns <- as.numeric(columns)
  data <- exprs(es)[rows, columns]

  ind <- which(is.na(data), arr.ind = T)
  if (nrow(ind) > 0) {
    data[ind] <- apply(data, 1, replacena, na.rm = T)[ind[,1]]
  }
  data <- t(data)

  pca <- prcomp(data)
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)

  xs <- sprintf("PC%s", seq_along(explained))
  xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)

  pData <- pData(es)[!(rownames(pData(es)) %in% setdiff(rownames(pData(es)), rownames(pca$x))),, drop=F]

  if (size != "") {
    pData[[size]] <- as.numeric(pData[[size]])
  }


  if (label == "id" || label == "") {
    label <- "names"
  }

#   pp <- ggplot(data=cbind(as.data.frame(pca$x), pData, sampleNames(es)))
#   if (size != "" && colour != "") {
#
#     aes <- aes_string(x=xs[n1],
#                       y=xs[n2], colour=colour, size=size)
#
#   } else if (colour != "") {
#     aes <- aes_string(x=xs[n1],
#                       y=xs[n2], colour=colour)
#   } else if (size != "") {
#     aes <- aes_string(x=xs[n1],
#                       y=xs[n2], size=size)
#   } else {
#     aes <- aes_string(x=xs[n1],
#                       y=xs[n2])
#   }

  pcadf <- as.data.frame(pca$x)
  names <- rownames(pcadf)
  gg <- plot_ly(data = cbind(pcadf, pData, names),
                type = "scatter",
                mode = "markers",
                x = ~eval(parse(text=xs[n1])),
                y = ~eval(parse(text=xs[n2])),
                color = if (colour != "") ~eval(parse(text=colour)) else 'rgba(0, 0, 0, .9)',
                marker = list(
                  size = if (size != "") ~eval(parse(text=size)) else 10),
                text = ~eval(parse(text=label))) %>%
        layout(xaxis = list(title = xlabs[n1], zeroline = F), yaxis = list(title = xlabs[n2], zeroline = F))

#   g <- pp + aes +
#     geom_point(aes) +
#     xlab(xlabs[n1]) + ylab(xlabs[n2])
#
#   #pg <- ggplotly(g)
#   if (label != "") {
#     message("i'm here 2s")
#
#     g <- g + aes + geom_text_repel(aes_string(label=label), size=3)
#   }
  #f <- tempfile(pattern="plot",tmpdir=getwd(),fileext=".svg")
  #ggsave(f, g)

  #print(capture.output(str(g)))
  return(tagList(gg))
}

