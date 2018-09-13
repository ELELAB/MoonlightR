#' FEA
#'
#' This function carries out the functional enrichment analysis (FEA)
#' @param BPname BPname biological process such as "proliferation of cells", "ALL" (default) if FEA should be carried out for all 101 biological processes
#' @param DEGsmatrix DEGsmatrix output from DEA such as dataDEGs"
#' @return matrix from FEA
#' @export
#' @examples
#' dataDEGs <- DPA(dataFilt = dataFilt,
#' dataType = "Gene expression")
#' dataFEA <- FEA(DEGsmatrix = dataDEGs)
FEA <- function (BPname = NULL, DEGsmatrix){

    DiseaseList <- get("DiseaseList")
    EAGenes <- get("EAGenes")

    if(is.null(BPname)){
        lf2 <- names(DiseaseList)
    }else{
        lf2 <- BPname
    }

    TableDiseasesNew <- NULL
    pb <- txtProgressBar(min = 0, max = length(DiseaseList), style = 3)

    for (k in 1:length(lf2)){
        setTxtProgressBar(pb, k)
        res <- as.data.frame(matrix(0,nrow = 1, ncol = 7, dimnames = list(1,c("Diseases.or.Functions.Annotation","p.Value",
                           "Moonlight.Z.score","commonNg","FunctionNg","Delta","Molecules") )))

        GeneList <- data.frame(PROBE_ID =  rownames(DEGsmatrix), logFC = DEGsmatrix$logFC, stringsAsFactors = FALSE)
        rownames(GeneList) <- GeneList$PROBE_ID

        selected_diseases <- as.data.frame(DiseaseList[[which(names(DiseaseList) == lf2[k])]])
        selected_diseases$ID <- selected_diseases$Genes.in.dataset

        res$commonNg <- length(intersect(GeneList$PROBE_ID, selected_diseases$ID))
        res$FunctionNg <- nrow(selected_diseases)
        res$Diseases.or.Functions.Annotation <- lf2[k]

        allgene <- unique(c(as.character(unique(EAGenes[, "ID"])), GeneList$PROBE_ID, selected_diseases$ID))

        seta <- allgene %in% GeneList$PROBE_ID
        setb <- allgene %in% selected_diseases$ID

        if( res$commonNg > 1){
        ft <- fisher.test(seta, setb)
        FisherpvalueTF <- ft$p.value
        res$p.Value <- FisherpvalueTF
        }
        else{ res$p.Value <- 1}

        GeneList <- GeneList[GeneList$PROBE_ID %in% selected_diseases$ID,]

        selected_diseases <- selected_diseases[selected_diseases$ID %in% GeneList[,"PROBE_ID"],]
        selected_diseases[,"Exp.Log.Ratio"] <- gsub(",",".",selected_diseases[,"Exp.Log.Ratio"])
        selected_diseases[,"Exp.Log.Ratio"] <- as.numeric(selected_diseases[,"Exp.Log.Ratio"])

        rownames(selected_diseases) <- selected_diseases$ID

        res$Molecules <- paste0(GeneList$PROBE_ID, collapse = ",")


        for( idx in 1: nrow(selected_diseases)){
            currTR <- selected_diseases$ID[idx]

            if (length(grep("Increases",selected_diseases[currTR,"Findings"]))==1){
                if(sign(GeneList[currTR,"logFC"]) > 0 ){
                    selected_diseases[currTR,"Prediction..based.on.expression.direction."]<-"Increased"
                }else if(sign(GeneList[currTR,"logFC"]) < 0 ){
                    selected_diseases[currTR,"Prediction..based.on.expression.direction."]<-"Decreased"
                }
            }
            if (length(grep("Decreases",selected_diseases[currTR,"Findings"]))==1){
                if(sign(GeneList[currTR,"logFC"]) < 0 ){
                    selected_diseases[currTR,"Prediction..based.on.expression.direction."]<-"Increased"
                }else if(sign(GeneList[currTR,"logFC"]) > 0 ){
                    selected_diseases[currTR,"Prediction..based.on.expression.direction."]<-"Decreased"
                }
            }
        }

        Correlation <- matrix(0,nrow(selected_diseases),1)
        selected_diseases <- cbind( selected_diseases, Correlation )

        if(length(grep("Decreases",selected_diseases$Findings)) != 0){

            selected_diseases[grep("Decreases",selected_diseases$Findings),"Findings"] <- -1
            selected_diseases[grep("Increases",selected_diseases$Findings),"Findings"] <- 1
            selected_diseases[grep("Affects",selected_diseases$Findings),"Findings"] <- 0
            selected_diseases[,"Findings"] <- as.numeric(selected_diseases[,"Findings"])

            PredictionIncreased <- which(sign(selected_diseases$Exp.Log.Ratio) == selected_diseases$Findings)
            PredictionDecreased <- which(sign(selected_diseases$Exp.Log.Ratio) != selected_diseases$Findings)
            PredictionAffected <- which(sign(selected_diseases$Findings) == 0)

            selected_diseases[PredictionIncreased,"Correlation"] <- 1
            selected_diseases[PredictionDecreased,"Correlation"] <- -1
            selected_diseases[PredictionAffected,"Correlation"] <- 0

            Zscore <- sum(selected_diseases$Correlation) / sqrt( length(PredictionIncreased) + length(PredictionDecreased) )
        }else{
            Zscore<-0
        }

        res$Moonlight.Z.score <- Zscore
        TableDiseasesNew <- rbind(TableDiseasesNew, res)
    }

    close(pb)

    TableDiseasesNew <- cbind(TableDiseasesNew, FDR = p.adjust(TableDiseasesNew$p.Value,method = "fdr") )

    return(TableDiseasesNew)
}


