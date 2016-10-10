#' @title plotNetworkHive: Hive network plot
#' @description
#' This function visualizes the GRN as a hive plot
#' @param dataGRN output GRN function
#' @param namesGenes list TSG and OCG to define axes
#' @param thres threshold of edges to be included
#' @import HiveR
#' @importFrom graphics plot.new
#' @importFrom grid gpar
#' @importFrom grDevices dev.cur
#' @export
#' @return no results Hive plot is executed
#' @examples
#' data(knownDriverGenes)
#' data(dataGRN)
#' plotNetworkHive(dataGRN = dataGRN, namesGenes = knownDriverGenes, thres = 0.55)
plotNetworkHive <- function(dataGRN, namesGenes, thres){

    names.genes.all <- intersect(as.character(unique(c(unlist(namesGenes), rownames(dataGRN[[1]])))),colnames(dataGRN[[1]]))
    tmp <- dataGRN[[1]][,(names.genes.all)]
    tmp[tmp<thres] <- 0

    genes.missing <- setdiff(names.genes.all, rownames(dataGRN[[1]]))
    tmp <- rbind(tmp, matrix(0, nrow=length(genes.missing), ncol=length(names.genes.all)))
    rownames(tmp)<- c(rownames(dataGRN[[1]]), genes.missing)

    tmp <- tmp[names.genes.all, names.genes.all]
    diag(tmp) <- 0

    myadj <- adj2HPD(M=tmp, axis.cols="lightgray")
    myadj$nodes$axis <- as.integer(rep(1, length(names.genes.all))+c(2*as.numeric(names.genes.all%in%namesGenes$TSG) +  as.numeric(names.genes.all%in%namesGenes$OCG)))
    n.axis <-  table(myadj$nodes$axis)
    names(n.axis) <- paste0("a.", names(n.axis))
    mycols <- c("darkgrey", "darkgreen","goldenrod")
    myadj$nodes$color <- mycols[myadj$nodes$axis]
    myadj$nodes$size <- 0.1
    for(i in 1:nrow(myadj$nodes)){
        myadj$nodes$radius[i] <- n.axis[paste0("a.",myadj$nodes$axis[i])]
        n.axis[paste0("a.",myadj$nodes$axis[i])] <- n.axis[paste0("a.",myadj$nodes$axis[i])]-1
    }

    ind.ocg <- which(rownames(tmp)[myadj$edges$id1] %in% namesGenes$OCG)
    myadj$edges$color[ind.ocg] <- "darkgreen"
    ind.ocg <- which(rownames(tmp)[myadj$edges$id2] %in% namesGenes$OCG)
    myadj$edges$color[ind.ocg] <- "darkgreen"

    ind.tsg <- which(rownames(tmp)[myadj$edges$id1] %in% namesGenes$TSG)
    myadj$edges$color[ind.tsg] <- "goldenrod"
    ind.tsg <- which(rownames(tmp)[myadj$edges$id2] %in% namesGenes$TSG)
    myadj$edges$color[ind.tsg] <- "goldenrod"

    pdf("networkHive.pdf")
    HiveR::plotHive(myadj, axLabs = c("remaining TFs", "OCG", "TSG"), bkgnd="white", anNode.gpar=gpar(fontsize = 10, col = "black", lwd = 0.5))
    if( (which = dev.cur()) != 1){graphics.off() }
}