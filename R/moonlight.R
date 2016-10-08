#' @title moonlight pipeline
#' @description moonlight is a tool for identification of cancer driver genes. 
#' This function wraps the different steps of the complete analysis workflow.
#' Providing different solutions:
#' \enumerate{
#' \item MoonlighR::FEA
#' \item MoonlighR::URA
#' \item MoonlighR::PIA
#' }
#' @param cancerType select cancer type for which analysis should be run. panCancer
#' 		for all available cancer types in TCGA. Defaults to panCancer
#' @param dataType dataType
#' @param directory directory 
#' @param BPname biological processes to use, if NULL: all processes will be used in analysis, RF for candidate; if not NULL the candidates for these processes will be determined (no learning)
#' @param cor.cut cor.cut Threshold
#' @param qnt.cut qnt.cut Threshold
#' @param Genelist Genelist
#' @param fdr.cut fdr.cut Threshold
#' @param logFC.cut logFC.cut Threshold
#' @param kNearest kNearest
#' @param nTF nTF
#' @param corThreshold corThreshold
#' @param nGenesPerm nGenesPerm
#' @param nBoot nBoot
#' @param DiffGenes DiffGenes
#' @param nSample nSample
#' @param thres.role thres.role
#' @param stage stage
#' @param subtype subtype
#' @param samples samples
#' @export
#' @return table with cancer driver genes TSG and OCG.
#' @examples
#' dataDEGs <- DPA(dataFilt = dataFilt, dataType = "Gene expression")
#' # to change with moonlight

moonlight <- function(cancerType="panCancer", dataType="Gene expression", 
                      directory =  "GDCdata", BPname = NULL,cor.cut = 0.6, 
                      qnt.cut = 0.25, Genelist= NULL, fdr.cut = 0.01, logFC.cut = 1,
                      corThreshold = 0.6, kNearest = 3, nGenesPerm = 10, DiffGenes = FALSE,
                      nBoot = 100, nTF = NULL, nSample=NULL,thres.role = 0, 
                      stage = NULL,subtype = 0, samples = NULL){

    if(length(cancerType) == 1 && cancerType == "panCancer"){
        cancerType <- sort(sapply(strsplit(grep("TCGA",GDCprojects,value=TRUE),"TCGA-"),"[",2))
    }

    res <- NULL

    for(cancer.i in cancerType){
        ### get TCGA data
        print("-----------------------------------------")
        print("Get TCGA data")
        print("-----------------------------------------")
        print(paste("cancer type:", cancer.i))

        dataFilt <- getDataTCGA(cancerType = cancer.i, dataType = dataType, 
                                  directory = directory, cor.cut = cor.cut, qnt.cut = qnt.cut, 
                                  nSample = nSample,stage = stage, 
                                  subtype = subtype, samples = samples)

        ### differential phenotype analysis
        print("-----------------------------------------")
        print("Differential phenotype analysis")
        print("-----------------------------------------")
        	
        dataDEGs <- DPA(dataType = dataType, dataFilt = dataFilt, fdr.cut = fdr.cut, 
        	                logFC.cut = logFC.cut)


        ### functional enrichment analysis -> carried out in URA, not necessary hear
        # print("-----------------------------------------")
        # print("Functional enrichment analysis")
        # print("-----------------------------------------")
        	
        # dataFEA <- FEA(BPname=BPname, DEGsmatrix = dataDEGs)

        ### gene regulatory network
        print("-----------------------------------------")
        print("Gene regulatory network")
        print("-----------------------------------------")

        #### parameter nTF for testing purposes
        if(is.null(nTF)){
            nTF <- nrow(dataDEGs)
        }

        if(is.null(Genelist)){
            Genelist <- rownames(dataDEGs)[1:nTF]
        }
        dataGRN <- GRN(TFs = Genelist, normCounts = dataFilt,
                       DEGsmatrix = dataDEGs,DiffGenes = FALSE,
                       nGenesPerm = nGenesPerm, kNearest = kNearest, nBoot = nBoot)

        ### upstream regulator analysis
        print("-----------------------------------------")
        print("Upstream regulator analysis")
        print("-----------------------------------------")

        dataURA <- URA(dataGRN = dataGRN, DEGsmatrix = dataDEGs, BPname = BPname)

        ### get TSG/OCG candidates using random forest
        print("-----------------------------------------")
        print("Get candidates")
        print("-----------------------------------------")
        listCandidates <- PRA(dataURA = dataURA, BPname = BPname, thres.role = thres.role)

        res.i <- list(dataDEGs = dataDEGs,
                dataURA = dataURA, 
                listCandidates = listCandidates)
        res <- c(res, list(res.i))
    }
    names(res) <- cancerType
    return(res)
}