 #' GSEA
#'
#' This function carries out the GSEA enrichment analysis.
#' @param DEGsmatrix DEGsmatrix output from DEA such as dataDEGs
#' @param top is the number of top BP to plot
#' @param plot if TRUE return a GSEA's plot
#' @importFrom grDevices dev.list
#' @importFrom grDevices graphics.off
#' @importFrom clusterProfiler bitr 
#' @importFrom DOSE gseDO 
#' @importFrom DOSE plot
#' @return return GSEA result
#' @export
#' @examples
#' dataDEGs <- DEGsmatrix
#' # dataFEA <- GSEA(DEGsmatrix = dataDEGs)
GSEA <- function (DEGsmatrix, top, plot = FALSE){

    dataDEGsnew <- cbind(mRNA = rownames(DEGsmatrix), DEGsmatrix)
  
  
  eg = as.data.frame(bitr(dataDEGsnew$mRNA,
                          fromType="SYMBOL",
                          toType="ENTREZID", 
                          OrgDb="org.Hs.eg.db"))
  eg <- eg[!duplicated(eg$SYMBOL),]
  dataDEGsFiltLevel <- dataDEGsnew
  dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
  
  dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
  eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]
  
  #all(eg$SYMBOL == dataDEGsFiltLevel$mRNA)
  dataDEGsFiltLevel$GeneID <- eg$ENTREZID
  
  dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
  genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
  names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
  
  genelistDEGs_sort <- sort(genelistDEGs,decreasing = TRUE)

  y <- gseDO(genelistDEGs_sort,
             nPerm         = 100, 
             minGSSize     = 120,
             pvalueCutoff  = 0.2, 
             pAdjustMethod = "BH",
             verbose       = FALSE)
  
  res <- as.matrix(summary(y))

  if (plot == TRUE){
  topID <- res[1,1]
 pdf("GSEAplot.pdf")
  DOSE::plot(y, geneSetID = topID)
  if (!(is.null(dev.list()["RStudioGD"]))){graphics.off()}
  }

  return(res)
  }


