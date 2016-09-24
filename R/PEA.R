#' PEA
#'
#' This function carries out the pathway enrichment analysis (FEA)
#' @param listMoonlight output Moonlight function
#' @importFrom clusterProfiler compareCluster
#' @importFrom clusterProfiler bitr 
#' @export
#' @return no return value, PEA result is plotted
PEA <- function(listMoonlight){

    list.ocg.id <- listMoonlight
    list.genes <- list()

    for ( i in 1:length(names(list.ocg.id))){

        listCur <- list.ocg.id[[i]]
        
        listCur<- listCur$dataDEGs
        DEGsmatrix<- listCur
        dataDEGsnew <- cbind(mRNA = rownames(DEGsmatrix), DEGsmatrix)
        eg = as.data.frame(bitr(dataDEGsnew$mRNA,
                                fromType="SYMBOL",
                                toType="ENTREZID", 
                                OrgDb="org.Hs.eg.db"))
        eg <- eg[!duplicated(eg$SYMBOL),]
       list.genes[[length(list.genes)+1]] <- eg$ENTREZID
         }
    
    names(list.genes)<- names(listMoonlight)
    # data(gcSample)
    res <- compareCluster(list.genes, fun="enrichPathway")
    plot(res)
    dev.off()

}
