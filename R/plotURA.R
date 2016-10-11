#' @title plotURA: Upstream regulatory analysis heatmap plot
#' @description
#' This function visualizes the URA in a heatmap
#' @param dataURA output URA function
#' @param plotNAME figure name
#' @importFrom gplots heatmap.2
#' @importFrom gplots redblue
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @return heatmap
#' @export
#' @examples 
#' data(dataURA)
#' dataDual <- PRA(dataURA = dataURA,
#' BPname = c("apoptosis","proliferation of cells"),
#' thres.role = 0)
#' plotURA(dataURA = dataURA[c(names(dataDual$TSG), names(dataDual$OCG)),],plotNAME = "URAplot")
plotURA<- function(dataURA, plotNAME = "URAplot"){    
    if(nrow(dataURA)>70){
        cexRow <- 0.1 + 1/(2*nrow(dataURA))*(log10(nrow(dataURA)))
    }else{
        cexRow <- 0.2 +  1/(10*log10(nrow(dataURA)))
    }

    pdf(file = paste0(plotNAME,".pdf"))
  
    par(oma=c(6,4,4,2))
    heatmap.2(dataURA, trace="none", col=rev(gplots::redblue(128)), Colv = TRUE, dendrogram = "row", 
    	 notecex=10/nrow(dataURA), cexCol =0.2 + 1/(3*log10(ncol(dataURA))), cexRow=cexRow,
    	 cellnote= signif(dataURA,4), notecol = "black",key=TRUE, keysize=2)   

    # if( (which = dev.cur()) != 1){ graphics.off() }
    # dev.off(which = dev.cur())
    dev.off()
    # dev.off()

}