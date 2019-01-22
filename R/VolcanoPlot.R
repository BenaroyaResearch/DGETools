#
# Volcano Plot, in progress
# Scott Presnell, SPresnell@benaroyaresearch.org
# Started spring 2018, Formalized July 13th, 2018
#
# Create a Volcano Plot
# run theme_set(theme_bw() + theme(text = element_text(face=fontFace, size=14), plot.title = element_text(face=fontFace, hjust = 0.5)))
# as part of your setup and repel font face will take from the current theme fontface
#
#' Create a Volcano plot
#'
#' Make a volcano plot, using cuts, ggrepel, optionally write out as image, optionally write out DEG list and all gene list to .csv
#'
#' @param fit fit ebays based fit after processing with limma
#' @param target string Name of the constrast target to plot
#' @param title string Title for the plot
#' @param numLabels integer Number of gene labels to show on the plot
#' @param write boolean Write DEG list (using contrast target string)
#' @param saveImage boolean Write out an image
#' @param ext string Image type/extension
#' @param labelSource string Name of the column for gene names
#' @param repelSize integer Size of labels
#' @param pointSize integet Size of points
#'
#' @return dataFrame Data Frame of differentially expressed gene information, in the form of topTable().
#'
#' @examples
#' Volcano(fit, "SPvsDP", "Single Positive vs Double Positive")
#'
#'
#' @export
#' @import limma gglot2 ggrepel
#'
Volcano <- function(fit,
                    target,
                    title,
                    numLabels=50,
                    write=FALSE,
                    saveImage=FALSE,
                    ext="png",
                    cutFC = 1.5,
                    cutPvalue = 0.05,
                    xlim=NULL,
                    labelSource="hgnc_symbol",  # mgi_symbol for mouse
                    repelFontFace = theme_get()$text$face,  # get the fontFace from the current theme by default.
                    repelSize=4,
                    pointSize=3) {


  # get the differentially expressed genes (actually get all genes, but sort them by p.value)
  tt <- topTable(fit, coef=target, sort.by="P", num=Inf)

  # Old way, use if gene list not provided to DEGList()
  #ttAnnot <- cbind("hgnc_symbol" = gene_key$hgnc_symbol[match(rownames(tt), gene_key$ensembl_gene_id)], tt)

  ttAnnot <- tt
  ttAnnot$lop = -log10(ttAnnot$adj.P.Val) # easier to keep in the dataframe
  ttAnnot$expr <- ifelse(ttAnnot$logFC >= log2(cutFC) & ttAnnot$adj.P.Val < cutPvalue, "up",
                         ifelse(ttAnnot$logFC < -log2(cutFC) & ttAnnot$adj.P.Val < cutPvalue, "down", "non"))

  # Only the differentially expressed ones for labeling, and possibly writing out.
  deg <- subset(ttAnnot, expr != "non")

  p <- ggplot(ttAnnot, aes(logFC, lop)) +
    geom_point(aes(color = expr), size=pointSize) +
    scale_colour_manual(name="Expression",
                        breaks = c("up", "down", "non"),
                        label  = c("Up regulated", "Down regulated", "Not significant"),
                        values = c("up"="red", "down"="blue", "non"="grey"))
  # Add an xlimit via xlim=xlim(-5,5), on the Volcano call.
  if (!is.null(xlim)) {
    p <- p + xlim
  }

  if (nrow(deg) > numLabels) {
    p <- p + geom_text_repel(data=deg[1:numLabels,], aes(logFC, lop, fontface=repelFontFace, label=deg[1:numLabels, labelSource]), size=repelSize)
  } else if (nrow(deg) > 0) {
    p <- p + geom_text_repel(data=deg, aes(logFC, lop, fontface=repelFontFace, label=deg[, labelSource]), size=repelSize)
  }

  p <- p +
    # do this outside
    # theme_bw() + theme(text = element_text(face=fontFace, size=14), plot.title = element_text(face=fontFace, hjust = 0.5)) +
    xlab(expression(log[2]~'Fold Change')) + ylab(expression(-log[10]~'p-value')) +
    ggtitle(title) +
    geom_hline(yintercept = -log10(cutPvalue), linetype = "dashed") +
    geom_vline(xintercept =  log2(cutFC),      linetype = "dashed") +
    geom_vline(xintercept = -log2(cutFC),      linetype = "dashed")

  print(p)

  if (saveImage == TRUE) {
    if (ext == "svg") {
      require(svglite)
    }
    # ggsave grabs the last plot if not specified
    ggsave(file=paste(target, ext, sep="."), width=8, height=10)
  }

  # pepare to write and/or return.
  if (nrow(deg) > 0) {
    deg$lop <- NULL
    deg$B <- NULL
  }
  #if (write == TRUE && nrow(deg) > 0) {
    # Row.names should be false if the Ensembl IDs are already added as a column possibly using DEGList(gene=) argument)
 #   write.csv(deg, paste(target, "deg.csv", sep="_"), quote=F)
 # } else if (write == TRUE) {
  if (write == TRUE) {
    ttAnnot$lop <- NULL
    ttAnnot$B <- NULL
    write.csv(ttAnnot, paste(target, "all.csv", sep="_"), quote=F)
  }

  # return only the differentially expressed genes.
  # invisible keeps the variable from being printed if the function call result is not assigned
  invisible(deg)
}
