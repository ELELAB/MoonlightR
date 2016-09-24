#' @title URA Upstream Regulator Analysis
#' @description
#' This function carries out the upstream regulator analysis
#' @param dataGRN output GNR function
#' @param DEGsmatrix output DPA function
#' @param BPname biological processes
#' @importFrom stats fisher.test
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
URA <- function(dataGRN, DEGsmatrix, BPname){
    if(is.null(BPname)){
        BPname <- names(DiseaseList)
    }
    #tRlist<- intersect(rownames(DEGsmatrix), rownames(dataGRN$miTFGenes))
    tRlist <- rownames(dataGRN$miTFGenes)
     # lf <- names(DiseaseList)

    TableDiseases <- matrix(0, nrow = length(tRlist), ncol = length(BPname), dimnames = list(tRlist, BPname))

    pb <- txtProgressBar(min = 0, max = nrow(TableDiseases), style = 3)

    for(j in 1:nrow(TableDiseases)) {
        currentTF <- as.character(rownames(TableDiseases)[j] )
        currentTF_regulon <- names(which(dataGRN$miTFGenes[currentTF,] > as.numeric(dataGRN$maxmi[currentTF])))
        currentTF_regulon <- as.matrix(currentTF_regulon)
        DEGsregulon <- intersect(rownames(DEGsmatrix), currentTF_regulon)
        if(length(DEGsregulon) > 2){
        tabFEA <- FEA(BPname = BPname, DEGsmatrix = DEGsmatrix[DEGsregulon,])
        TableDiseases[j,] <- tabFEA$Activation.z.score
        }
        else { TableDiseases[j,] <- 0}
    }

    close(pb)

    return(TableDiseases)

}