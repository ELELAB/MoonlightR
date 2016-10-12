#' @title Generate network 
#' @description This function carries out the gene regulatory network inference using parmigene
#' @param TFs a vector of genes.
#' @param DEGsmatrix DEGsmatrix output from DEA such as dataDEGs
#' @param DiffGenes if TRUE consider only diff.expr genes in GRN
#' @param normCounts is a matrix of gene expression with genes in rows and samples in columns.
#' @param kNearest the number of nearest neighbors to consider to estimate the mutual information.
#' @param seed set to get same result
#' Must be less than the number of columns of normCounts.
#' @param nGenesPerm nGenesPerm
#' @param nBoot nBoot
#' @importFrom parmigene knnmi.cross
#' @export
#' @return an adjacent matrix
#' @examples
#' dataDEGs <- DEGsmatrix
#' dataGRN <- GRN(TFs = rownames(dataDEGs)[1:100], 
#' DEGsmatrix = dataDEGs,
#' DiffGenes = TRUE,
#' normCounts = dataFilt)
GRN <- function(TFs, DEGsmatrix, DiffGenes = FALSE, normCounts, kNearest = 3, nGenesPerm = 10, nBoot = 10, seed=12345) {
    set.seed(seed)
    normCountsA <- normCounts
    normCountsB <- normCounts

    if(DiffGenes==TRUE){
      commonGenes <- intersect(rownames(DEGsmatrix), rownames(normCountsB) )
      normCountsB <- normCountsB[commonGenes,]
    }else{
        normCountsB <- normCountsA
    }


    MRcandidates <- intersect(rownames(normCountsA),TFs) 


    # Mutual information between TF and genes
    sampleNames <- colnames(normCounts)
    geneNames <- rownames(normCounts)

    # messageMI_TFgenes <- paste("Estimation of MI among [", length(MRcandidates), " TRs and ", nrow(normCounts), " genes].....", sep = "")
    # timeEstimatedMI_TFgenes1 <- length(MRcandidates)*nrow(normCounts)/1000
    # timeEstimatedMI_TFgenes <- format(timeEstimatedMI_TFgenes1*ncol(normCounts)/17000, digits = 2)
    # messageEstimation <- print(paste("I Need about ", timeEstimatedMI_TFgenes, "seconds for this MI estimation. [Processing 17000k elements /s]  "))

    # system.time(
    miTFGenes <- knnmi.cross(normCountsA[MRcandidates, ], normCountsB, k = kNearest)
    # )


    # threshold with bootstrap   
    tfListCancer <- TFs
    tfListCancer <- intersect(tfListCancer,rownames(normCountsA))

    maxmi<-rep(0,length(tfListCancer))

    Cancer_null_distr<-matrix(0,length(tfListCancer),nBoot)
    rownames(Cancer_null_distr)<-tfListCancer

    for (i in 1: nBoot){
    # cat(paste( (nBoot-i),".",sep=""))
        SampleS <- sample(1:ncol(normCountsA))
        g <- sample(1:nrow(normCountsA), nGenesPerm)
        # if(i == 1) system.time(mi <- knnmi.cross(normCounts[tfListCancer, ], normCounts[g, SampleS], k = kNum)) else
        mi <- knnmi.cross(normCountsA[tfListCancer, ], normCountsA[g, SampleS], k = kNearest)

        maxmiCurr <- apply(mi,1, max)
        Cancer_null_distr[,i] <- maxmiCurr
        index <- maxmi < maxmiCurr
        maxmi[index]<- maxmiCurr[index]
    }

    names(maxmi) <- rownames(Cancer_null_distr)


    return(list(miTFGenes = miTFGenes, maxmi = maxmi))
  
}