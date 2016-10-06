#' @title URA Upstream Regulator Analysis
#' @description
#' This function carries out the upstream regulator analysis
#' @param dataGRN output GNR function
#' @param DEGsmatrix output DPA function
#' @param BPname biological processes
#' @param nCores number of cores to use
#' @importFrom stats fisher.test
#' @import doParallel
#' @import foreach
#' @export
#' @return an adjacent matrix
#' @examples
#' dataDEGs <- DEGsmatrix
#' dataGRN <- GRN(TFs = rownames(dataDEGs)[1:100], 
#' DEGsmatrix = dataDEGs,
#' DiffGenes = TRUE,
#' normCounts = dataFilt)
#' dataURA <-URA(dataGRN = dataGRN,
#' DEGsmatrix = dataDEGs, 
#' BPname = c("apoptosis",
#' "proliferation of cells"))
URA <- function(dataGRN, DEGsmatrix, BPname, nCores = 1){
doParallel::registerDoParallel(cores = nCores)
 
  
    if(is.null(BPname)){
        BPname <- names(DiseaseList)
    }
    #tRlist<- intersect(rownames(DEGsmatrix), rownames(dataGRN$miTFGenes))
    tRlist <- rownames(dataGRN$miTFGenes)
     # lf <- names(DiseaseList)

    pb <- txtProgressBar(min = 0, max = length(tRlist), style = 3)

    TableDiseases <- foreach(j = 1:length(tRlist), .combine = "rbind", .packages="foreach") %dopar% {
      setTxtProgressBar(pb, j)
        currentTF <- as.character(tRlist[j] )
        currentTF_regulon <- names(which(dataGRN$miTFGenes[currentTF,] > as.numeric(dataGRN$maxmi[currentTF])))
        currentTF_regulon <- as.matrix(currentTF_regulon)
        DEGsregulon <- intersect(rownames(DEGsmatrix), currentTF_regulon)
        if(length(DEGsregulon) > 2){
        tabFEA <- FEA(BPname = BPname, DEGsmatrix = DEGsmatrix[DEGsregulon,])
        return(tabFEA$Activation.z.score)
        }else{
            return(rep(0, length(BPname)))
        }
    }
    dimnames(TableDiseases) <- list(tRlist, BPname)
 
    close(pb)
    return(TableDiseases)
}