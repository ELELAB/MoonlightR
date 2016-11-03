#' @title plotFEA
#' @description 
#' This function visualize the functional enrichment analysis (FEA)'s barplot
#' @param dataFEA dataFEA 
#' @param topBP topBP 
#' @param additionalFilename additionalFilename 
#' @param height Figure height
#' @param width Figure width
#' @param offsetValue offsetValue 
#' @param angle angle 
#' @param xleg xleg 
#' @param yleg yleg 
#' @param bar.colors Colors of the bars.
#' @param minY minY 
#' @param maxY maxY 
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @return no return value, FEA result is plotted
#' @examples
#' dataFEA <- FEA(DEGsmatrix = DEGsmatrix)
#' plotFEA(dataFEA = dataFEA, additionalFilename = "_example",height = 20,width = 10)
plotFEA <- function(dataFEA, 
                    topBP = 10, 
                    additionalFilename = NULL, 
                    height, 
                    width,
                    bar.colors = c("darkgreen", "gold", "navy"),
                    offsetValue = 5, 
                    angle = 90, 
                    xleg = 35, 
                    yleg =5 , 
                    minY = -5, 
                    maxY =10){
  
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
  
  toPlot[1, which (toPlot[1, ] > 10)] <- 10 # check on y limits
  
  #toPlot <-as.numeric(toPlot)
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
  }
  colnames(toPlot) <-  paste(sapply(tmp$Diseases.or.Functions.Annotation, simpleCap), " (n=", tmp$commonNg, ")")
  rownames(toPlot) <- c("-logFDR/10","Activation Z-score Increased","Activation Z-score Decreased")
  
  toPlot <- melt(toPlot)
  colnames(toPlot) <- c("score","Biofunction","value")
  p <- ggplot(toPlot, aes(x=Biofunction, y=value, fill=factor(score))) +
    geom_bar(stat="identity", position="dodge", colour="white") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("FEA - Enriched BioFunctions") + theme(plot.title = element_text(lineheight=.8, face="bold")) +  xlab("") +  ylab("-logFDR/10 and Activation z-score") + labs(fill="") +
    scale_fill_manual(values=bar.colors)  +  facet_grid(.~Biofunction,scales="free_x") + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.y=element_text(margin=margin(0,20,0,10)),
          axis.text.x = element_text(size = 12),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          legend.key = element_rect(colour = 'white'),
          #axis.text.y = element_text(angle = 90, hjust = 1),
          legend.position="top",
          axis.ticks.y = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          panel.margin.x=unit(1.5, "lines"),
          strip.text.y = element_blank()
    ) +  scale_y_continuous(breaks=c(seq(minY,maxY,1)),expand = c(0, 0.6)) +
    coord_cartesian(ylim = c(min(toPlot$value),max(toPlot$value))) +
    geom_hline(aes(yintercept = 0),
               colour = "black")
  

  if(!is.null(additionalFilename)){
    ggsave(plot = p, file = paste0("plotFEA",additionalFilename,".pdf"))
  } else {
    plot(p)
  }
  
}