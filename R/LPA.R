#' LPA
#'
#' This function carries out the literature phenotype analysis (LPA)
#' @param dataDEGs is output from DEA
#' @param BP is biological process
#' @param BPlist is list of genes annotated in BP
#' @import RISmed 
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#' @return table with number of pubmed that affects, increase or decrase genes annotated in BP
#' @examples
#' data(DEGsmatrix)
#' BPselected <- c("apoptosis")
#' BPannotations <- DiseaseList[[match(BPselected, names(DiseaseList))]]$ID
#' dataLPA <- LPA(dataDEGs = DEGsmatrix[1:5,],
#'                  BP =  BPselected,
#'                  BPlist = BPannotations)
LPA <- function (dataDEGs, BP, BPlist) {
    
  BPgenesDEGs <- intersect(BPlist, rownames(dataDEGs))
  dataDEGsBP <- dataDEGs[BPgenesDEGs,]

  DiseaseMN <- matrix(0, nrow(dataDEGsBP),7)
  colnames(DiseaseMN) <- c("ID",
                               "Genes.in.dataset",
                               "Prediction..based.on.expression.direction.",
                               "Exp.Log.Ratio",
                               "Findings",
                               "PubmedDecreases",
                               "PubmedIncreases")
  DiseaseMN <- as.data.frame(DiseaseMN)
  DiseaseMN$Prediction..based.on.expression.direction. <- rep("Decreased", nrow(dataDEGsBP))
  
  DiseaseMN$ID <- BPgenesDEGs
  DiseaseMN$Genes.in.dataset <- BPgenesDEGs
  DiseaseMN$Exp.Log.Ratio <- dataDEGs[BPgenesDEGs,"logFC"]
  rownames(DiseaseMN) <- DiseaseMN$ID

  pb <- txtProgressBar(min = 0, max = nrow(DiseaseMN), style = 3)
  
  for ( i in 1: nrow(DiseaseMN)){
  
    setTxtProgressBar(pb, i)
    curG <- DiseaseMN$ID[i]
    
  keypubmed <- paste(curG,BP,"decreases")
  res <- EUtilsSummary(keypubmed, type="esearch", db="pubmed")
  fetch <- EUtilsGet(res,type="efetch", db="pubmed")
  RecordsDec <- length(fetch@PMID)
  DiseaseMN[curG,"PubmedDecreases"] <- RecordsDec
  
  keypubmed <- paste(curG,BP,"increases")
  res <- EUtilsSummary(keypubmed, type="esearch", db="pubmed")
  fetch <- EUtilsGet(res,type="efetch", db="pubmed")
  RecordsInc <- length(fetch@PMID)
  
  DiseaseMN[curG,"PubmedIncreases"] <- RecordsInc
   
  if(DiseaseMN[curG,"PubmedDecreases"] == DiseaseMN[curG,"PubmedIncreases"]){
    DiseaseMN[curG,"Findings"] <- paste0("Affects (", RecordsInc, ")")
  }
  
  if(DiseaseMN[curG,"PubmedDecreases"] < DiseaseMN[curG,"PubmedIncreases"]){
    DiseaseMN[curG,"Findings"] <- paste0("Increases (", RecordsInc, ")")
  }
  
  if(DiseaseMN[curG,"PubmedDecreases"] > DiseaseMN[curG,"PubmedIncreases"]){
    DiseaseMN[curG,"Findings"] <- paste0("Decreases (", RecordsDec, ")")
  }
  
  }
  close(pb)
  
  return(DiseaseMN)
}
