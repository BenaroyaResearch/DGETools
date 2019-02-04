#
# Volcano Plot, in progress
# Scott Presnell, SPresnell@benaroyaresearch.org
# Started spring 2018, Formalized July 13th, 2018
#
#
#if(getRversion() >= "3.1.0") utils::globalVariables(c("expr", logFC", "lop"))
# Create a Volcano Plot
# run theme_set(theme_bw() + theme(text = element_text(face=fontFace, size=14), plot.title = element_text(face=fontFace, hjust = 0.5)))
# as part of your setup and repel font face will take from the current theme fontface
#
#' Create a Volcano plot
#'
#' Make a volcano plot, using cuts, ggrepel, optionally write out as image, optionally write out DEG list and all gene list to .csv
#'
#' @param fit A fitted linear model object (MArrayLM)
#' @param target Name of the constrast target to plot (string)
#' @param title Title for the plot (string)
#' @param writeTopTable Write gene list to CSV file (see \code{target} string, boolean)
#' @param write deprecated in favor of \code{witeTopTable}
#' @param saveImage Write out an image (boolean)
#' @param extImage Image type/extension (string)
#' @param ext deprecated in favor of \code{extImage}
#' @param cutFC Fold Change cutoff for significance (floating point)
#' @param cutPvalue p-valut cutoff for significance (floating point)
#' @param numLabels  Number of labels to print (integer)
#' @param labelList List of label names to render, only those (vector)
#' @param labelSize Size of text in repel for labels from labelList - default is the same as the repelSize (integer)
#' @param labelSource Name of the column for gene names (string)
#' @param colors vector of up/right point color, down/left point color, and not significant color (string)
#' @param repelFontFace Name of fontface for labels (string)
#' @param repelSize Size of text repel labels (integer)
#' @param pointSize Size of geom_points (integer)
#' @param ...  Additional function calls from the ggplot2 package to be added to ggplot, such as \code{xlim=xlim(-5,5)}, \code{ylim=ylim(0,9)}
#'
#' @return A dataframe of differentially expressed gene information, in the form of topTable().
#'
#' @author Scott R Presnell, \email{SPresnell@@benaroyaresearch.org}
#' @seealso Associated or similar limma functions \code{\link[limma]{topTable}},  \code{\link[limma]{volcanoplot}}
#'
#' @examples
#' \donttest{
#' Volcano(fit, "SPvsDP", "Single Positive vs Double Positive")
#'}
#'
#'
#' @export
#' @import limma ggplot2 ggrepel
#' @importFrom utils write.csv
#'
Volcano <- function(fit,
                    target,
                    title,
                    numLabels=50,
                    write=FALSE,
                    writeTopTable=FALSE,
                    saveImage=FALSE,
                    ext="png",
                    extImage="png",
                    cutFC = 1.5,
                    cutPvalue = 0.05,
                    labelList=NULL,
                    labelSize = 4,
                    labelSource="hgnc_symbol",  # mgi_symbol for mouse
                    colors=c("red", "blue", "grey"), # right, left, not significant
                    repelFontFace = theme_get()$text$face,  # get the fontFace from the current theme by default.
                    repelSize=4,
                    pointSize=3,
                    ...) {


  if (!missing(write)) {
    warning("argument 'write' is deprecated, please use writeTopTable instead.", call. = FALSE)
    writeTopTable <- write
  }

  if (!missing(ext)) {
    warning("argument 'ext' is deprecated, please use extImage instead.", call. = FALSE)
    extImage <- ext
  }

  # copy of the repelSize to the labelSize if labelSize is not set.
  if (missing(labelSize)) {
    labelSize <- repelSize
  }

  # get all genes, sorted by p.value
  topResults <- topTable(fit, coef=target, sort.by="P", number=Inf)
  topResults$lop <- -log10(topResults$adj.P.Val) # easier to keep in the dataframe
  topResults$expr <- ifelse(topResults$logFC >= log2(cutFC) & topResults$adj.P.Val < cutPvalue, "up",
                         ifelse(topResults$logFC < -log2(cutFC) & topResults$adj.P.Val < cutPvalue, "down", "non"))
  topResults$labels <- as.character(topResults[[labelSource]])

  # Old way uses global variable,  used when the gene list was not provided to DEGList()
  #topResults$labels <- gene_key$hgnc_symbol[match(rownames(topResults), gene_key$ensembl_gene_id)], topResults)

  # Only the differentially expressed ones for labeling.
  degFrame <- subset(topResults, expr != "non")

  # If we want a limited set of labels, then find the indicies for those labels, and copy them onto a background of empty strings
  # (if we use NA, geomp_text_repel complains, and it doesn't calculate the "bumps" correctly)
  if (!is.null(labelList)) {
    newLabels <- rep("", length(topResults$labels)) # label text strings
    newSize <- rep(repelSize, length(topResults$labels)) # label text size
    indices <-  which(topResults$labels %in% labelList) #indices of desired labels
    nDEGCount = nrow(degFrame)
    maxLabels = nDEGCount

    #determine the number of labels to render.
    if (nDEGCount > numLabels) {
      maxLabels = numLabels
    }

    # combine top numLabels and the requested labels
    if (maxLabels > 0) {
      newLabels[c(1:maxLabels, indices)] <- topResults$labels[c(1:maxLabels, indices)]
    # ... or just the requested labels
    } else {
      newLabels[indices] <- topResults$labels[indices]
    }
    newSize[indices] <- labelSize
    topResults$labels <- newLabels
  }


  myPlot <- ggplot(data=topResults, aes_(x=quote(logFC), y=quote(lop))) +
    geom_point(aes_(color=quote(expr)), size=pointSize) +
    scale_colour_manual(name="Expression",
                        breaks = c("up", "down", "non"),
                        label  = c("Up regulated", "Down regulated", "Not significant"),
                        values = c("up"=colors[1], "down"=colors[2], "non"=colors[3]))

  # grab the variable argument list from the elipsis in the function definition.
  variableArguments <- list(...)

  # Loop over it to add elements to the plot - assumes variable arguments are only to be fed to ggplot()
  # If variableArguments is an empty list, loop doesn't happen, no error.
  for (arg in variableArguments) {
    myPlot <- myPlot + arg
  }

  if (!is.null(labelList)) {
    myPlot <- myPlot + geom_text_repel(data=topResults, aes(fontface=repelFontFace, label=labels), size=newSize)
  } else if (nrow(degFrame) > numLabels) {
    myPlot <- myPlot + geom_text_repel(data=degFrame[1:numLabels,], aes(fontface=repelFontFace, label=degFrame[1:numLabels, "labels"]), size=repelSize)
  } else if (nrow(degFrame) > 0) {
    myPlot <- myPlot + geom_text_repel(data=degFrame, aes(fontface=repelFontFace, label=labels), size=repelSize)
  }
    # fall through - no differentially expressed genes, hence no labels.

  myPlot <- myPlot +
    # do this outside the function
    # theme_set(theme_bw() + theme(text = element_text(face="bold", size=14), plot.title = element_text(face="bold", hjust = 0.5)))
    xlab(expression(log[2]~'Fold Change')) + ylab(expression(-log[10]~'p-value')) +
    ggtitle(title) +
    geom_hline(yintercept = -log10(cutPvalue), linetype = "dashed") +
    geom_vline(xintercept =  log2(cutFC),      linetype = "dashed") +
    geom_vline(xintercept = -log2(cutFC),      linetype = "dashed")

  print(myPlot)

  if (saveImage == TRUE) {
#    if (extImage == "svg") {
#    }
    # ggsave grabs the last plot if not specified
    ggsave(filename=paste(target, extImage, sep="."), width=8, height=10, dpi=320)
  }

  # write out the top table results to a CSV file with all the results, remove redunant information
  if (writeTopTable == TRUE) {
    topResults$lop <- NULL
    topResults$B <- NULL
    topResults$labels <- NULL
    write.csv(topResults, paste(target, "all.csv", sep="_"), quote=F)
  }

  # return only the differentially expressed genes.
  # invisible keeps the variable from being printed if the function call result is not assigned
  # pepare to return.
  if (nrow(degFrame) > 0) {
    degFrame$lop <- NULL
    degFrame$labels <- NULL
  }
  invisible(degFrame)
}
