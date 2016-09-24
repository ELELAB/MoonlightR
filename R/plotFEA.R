#' @title plotFEA
#' @description 
#' This function visualize the functional enrichment analysis (FEA)'s barplot
#' @param dataFEA dataFEA 
#' @param topBP topBP 
#' @param plotNAME plotNAME 
#' @param height Figure height
#' @param width Figure width
#' @param offsetValue offsetValue 
#' @param angle angle 
#' @param xleg xleg 
#' @param yleg yleg 
#' @param minY minY 
#' @param maxY maxY 
#' @importFrom graphics barplot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics text
#' @export
#' @return no return value, FEA result is plotted
#' @examples
#' dataFEA <- FEA(DEGsmatrix = DEGsmatrix)
#' plotFEA(dataFEA = dataFEA,plotNAME = "FEAplot",height = 20,width = 10)
plotFEA <- function(dataFEA, topBP = 10,plotNAME = "test",height,width,
                               offsetValue=5,angle=90,xleg=35,yleg=5,minY=-5,maxY =10){
                               #offsetValue,angle,xleg=22,yleg=10,minY=-3,maxY=10){
  
 # mycols <- c("#66C2A5", "#FC8D62", "#8DA0CB")
    mycols <- c("#8DD3C7", "#FFFFB3", "#BEBADA")
    pdf(file = paste0(plotNAME,".pdf"))
    par(mar = c(12,5,5,1))
    dataFEA <- dataFEA[order(abs(dataFEA$Activation.z.score),decreasing =TRUE),]

    tmp <- dataFEA[1:topBP,]

    # tmp <- chosenFile[toupper(chosenFile[,"Filter"]) == "YES",]

    tmp <- as.data.frame(tmp)
    tmp$p.Value <- as.numeric(tmp$p.Value)
    tmp$Activation.z.score <- as.numeric(tmp$Activation.z.score)
    tmp$commonNg <- as.numeric(tmp$commonNg)


    tmp$p.Value <- -log2(tmp$p.Value)/10


    toPlot <- matrix(0, nrow = 3, ncol = nrow(tmp))
    #toPlot <- as.data.frame(toPlot) 
    toPlot[1, ] <- tmp$p.Value
    toPlot <- toPlot[, order(toPlot[1, ], decreasing = TRUE)]
    toPlot[2, as.numeric(tmp$Activation.z.score) > 0] <- as.numeric(tmp$Activation.z.score[as.numeric(tmp$Activation.z.score) > 0])
    toPlot[3, as.numeric(tmp$Activation.z.score) < 0] <- as.numeric(tmp$Activation.z.score[as.numeric(tmp$Activation.z.score) < 0])


    toPlot[1, which (toPlot[1, ] > 10)] <-10 # check on y limits

    #toPlot <-as.numeric(toPlot)
    toPlot <- as.matrix(toPlot)
    xAxis <- barplot(toPlot, beside = TRUE, col = mycols, ylab = "-logFDR/10 and Activation z-score",
                   names = NULL, main = paste(plotNAME, " - Enriched BioFunctions",sep=" "), ylim = c(minY, maxY))
    legend("topright", bty = "n", xleg, yleg, legend = c("-logFDR/10","Activation Z-score Increased","Activation Z-score Decreased"), 
         fill = mycols, cex=0.5)#, text.col = c("blue","red", "green"), pch = 15)

    text(xAxis[2, ], par("usr")[3]+offsetValue,  srt = angle, adj = 1, xpd = TRUE,
       labels = paste(tmp$Diseases.or.Functions.Annotation, " (n=", tmp$commonNg, ")", sep = ""), cex = 1)


    dev.off()
}