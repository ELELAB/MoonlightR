#' DPA
#'
#' This function carries out the differential phenotypes analysis
#' @param dataConsortium is TCGA or GEO, default TCGA
#' @param dataFilt  obtained from get.data.TCGA
#' @param dataType  selected
#' @param fdr.cut is a threshold to filter DEGs according their p-value corrected
#' @param logFC.cut is a threshold to filter DEGs according their logFC
#' @param diffmean.cut diffmean.cut for DMR
#' @param samplesType samplesType
#' @param colDescription colDescription
#' @param gset gset
#' @param gsetFile gsetFile
#' @importFrom TCGAbiolinks TCGAanalyze_DEA
#' @importFrom TCGAbiolinks TCGAanalyze_DMR
#' @importFrom TCGAbiolinks TCGAquery_SampleTypes
#' @import SummarizedExperiment
#' @import Biobase 
#' @importFrom stats quantile
#' @importFrom stats model.matrix
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @return result matrix from differential phenotype analysis
#' @export
#' @examples 
#' dataDEGs <- DPA(dataFilt = dataFilt, dataType = "Gene expression")

DPA <- function(dataType, dataFilt, dataConsortium="TCGA", fdr.cut = 0.01, logFC.cut = 1, 
                diffmean.cut = 0.25, samplesType, colDescription, gset, 
                gsetFile = "gsetFile.RData"){
  
    if (dataConsortium == "GEO"){
        tabGset <- as.data.frame(subset(pData(gset), select=c("geo_accession",colDescription,"type")))
        colnames(tabGset) <- c("SampleID", "Disease","Group")
        tabGset$SampleID <- as.character(tabGset$SampleID)
        tabGset$Disease <- as.character(tabGset$Disease)
        tabGset$Group <- rep("X", length(tabGset$Group ))

        tabGset[grep (tolower(samplesType[1]), tolower(tabGset$Disease)), "Group"]<-0
        tabGset[grep (tolower(samplesType[2]), tolower(tabGset$Disease)), "Group"]<-1
        gsms <- paste(tabGset$Group,collapse = "")

        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        sel <- which(sml != "X")
        sml <- sml[sel]
        gset <- gset[ ,sel]

        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
        LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }

        sml <- paste("G", sml, sep="")    
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2,adjust.method ="fdr", sort.by="B", number=50000)

        #save(gset, tT, tabGset, file = gsetFile)

        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

        Cond1num <- sum(tabGset$Group == 0)
        Cond2num <- sum(tabGset$Group == 1)

        message("----------------------- DEA -------------------------------")
        message(message1 <- paste("there are Cond1 type", samplesType[1], 
                                "in ", Cond1num, "samples"))
        message(message2 <- paste("there are Cond2 type", samplesType[2], 
                                "in ", Cond2num, "samples"))
        message(message3 <- paste("there are ", length( unique(tT$Gene.symbol)), "features as miRNA or genes "))


        tT_filt<- tT[tT$adj.P.Val < fdr.cut,]
        tT_filt_FC <- tT_filt[abs(tT_filt$logFC) >= logFC.cut,]

        tT_filt_FC <- tT_filt_FC[order(abs(tT_filt_FC$logFC), decreasing = TRUE),]
        tT_filt_FC <- tT_filt_FC[!duplicated(tT_filt_FC$Gene.symbol),]
        tT_filt_FC <- tT_filt_FC[tT_filt_FC$Gene.symbol!="",]
        tT_filt_FC <- tT_filt_FC[!is.na(tT_filt_FC$Gene.symbol),]


        rownames(tT_filt_FC) <- tT_filt_FC$Gene.symbol

        dataDEGs<- tT_filt_FC
        return(dataDEGs)  
  
    }else if(dataConsortium == "TCGA"){
        if(dataType == "Gene expression"){

            dataSmTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "TP")
            dataSmNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = "NT")

            dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                                  mat2 = dataFilt[,dataSmTP],
                                  Cond1type = "Normal",
                                  Cond2type = "Tumor",
                                  fdr.cut = fdr.cut ,
                                  logFC.cut = logFC.cut,
                                  method = "glmLRT")
            return(dataDEGs)
        }else if(dataType == "Methylation"){
            dataFilt <- subset(dataFilt,subset = (rowSums(is.na(assay(dataFilt))) == 0))
        
            system.time (cancer.met <- TCGAanalyze_DMR(dataFilt, groupCol = "shortLetterCode",
                                                 group1 = "TP",
                                                 group2= "NT",
                                                 p.cut = fdr.cut,
                                                 diffmean.cut = 0.25,
                                                 legend = "State",
                                                 plot.filename = "test.png"))
      
            values(dataDEGs)$status.NT.TP
            sig.met <- dataDEGs[values(dataDEGs)$status.NT.TP %in% c("Hypermethylated","Hypomethylated"),]
        
            return(cancer.met)
        }
    }
  
}