#' getDataTCGA
#'
#' This function retrieves and prepares TCGA data
#' @param cancerType select cancer type for which analysis should be run. panCancer
#' 		for all available cancer types in TCGA. Defaults to panCancer
#' @param dataType is dataType such as gene expression, cnv, methylation etc.
#' @param directory	Directory/Folder where the data was downloaded. Default: GDCdata
#' @param cor.cut cor.cut 
#' @param qnt.cut qnt.cut 
#' @param nSample nSample 
#' @param stage stage 
#' @param subtype subtype
#' @param samples samples
#' @param seed set to get same result
#' @importFrom TCGAbiolinks GDCquery
#' @importFrom TCGAbiolinks GDCquery_clinic
#' @importFrom TCGAbiolinks TCGAquery_SampleTypes
#' @importFrom TCGAbiolinks TCGAquery_subtype
#' @importFrom TCGAbiolinks GDCdownload
#' @importFrom TCGAbiolinks GDCprepare
#' @importFrom TCGAbiolinks TCGAanalyze_Preprocessing
#' @importFrom TCGAbiolinks TCGAanalyze_Normalization
#' @importFrom TCGAbiolinks TCGAanalyze_Filtering
#' @importFrom utils as.roman
#' @export
#' @return returns filtered TCGA data
#' @examples
#' \dontrun{
#' dataFilt <- getDataTCGA(cancerType = "LUAD", 
#' dataType = "Gene expression", directory = "data", nSample = 4)
#' }

getDataTCGA <- function(cancerType, dataType, directory, 
                          cor.cut = 0.6, qnt.cut = 0.25, 
                          nSample, stage = "ALL", subtype = 0, 
                          samples = NULL, seed = 12345){

  DiseaseList <- get("DiseaseList")
  GDCprojects <- get("GDCprojects")
  geneInfo <- get("geneInfo")
  
    set.seed(seed)
    CancerProject <- paste0("TCGA-",cancerType)
    DataDirectory <- paste0(directory,"GDC_",gsub("-","_",CancerProject))

    if(dataType == "Gene expression"){

        FileNameData <- paste0(DataDirectory, "_stage_", stage,
                               "_subtype_", subtype,
                               "_","IlluminaHiSeq",".rda")
        if(!file.exists(FileNameData)){
            query <- GDCquery(project = CancerProject, 
                              data.category = "Gene expression",
                              data.type = "Gene expression quantification",
                              platform = "Illumina HiSeq", 
                              file.type = "results",
                              legacy = TRUE)

            samplesDown <- query$results[[1]]$cases

            dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                              typesample = "TP")

            dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                              typesample = "NT")
            
            if (is.numeric(stage)==TRUE){
                dataClin <- GDCquery_clinic(project = CancerProject,type = "clinical_patient")
                curStage <- paste0("Stage ", as.roman(stage))
                dataClin$tumor_stage <- toupper(dataClin$tumor_stage)
                dataClin$tumor_stage <- gsub("[ABCDEFGH]","",dataClin$tumor_stage)
                dataClin$tumor_stage <- gsub("ST","Stage",dataClin$tumor_stage)

                dataStg <- dataClin[dataClin$tumor_stage %in% curStage,]
                message(paste(curStage, "with", nrow(dataStg), "samples"))

                dataStgC <- dataSmTP[substr(dataSmTP,1,12) %in% dataStg$bcr_patient_barcode]
                dataSmTP <- dataStgC
            }
            
            if (is.character(subtype)==TRUE){
              dataSmSubt <- dataSmTP[substr(dataSmTP,1,12) %in% samples]
              dataSmTP <- dataSmSubt
            }
            
            message <- paste0(CancerProject," ... with ", length(dataSmTP), " TP, ", length(dataSmNT), " NT")
            print(message)

            if(is.null(nSample)){
                queryDown <- GDCquery(project = CancerProject, 
                                  data.category = "Gene expression",
                                  data.type = "Gene expression quantification",
                                  platform = "Illumina HiSeq", 
                                  file.type = "results",
                                  barcode = c(dataSmTP, dataSmNT),
                                  # barcode = c(sample(dataSmTP,nSample), sample(dataSmNT,nSample)),
                                  legacy = TRUE)
            }else{
                queryDown <- GDCquery(project = CancerProject, 
                                  data.category = "Gene expression",
                                  data.type = "Gene expression quantification",
                                  platform = "Illumina HiSeq", 
                                  file.type = "results",
                                  # barcode = c(dataSmTP, dataSmNT),
                                  barcode = c(sample(dataSmTP,nSample), sample(dataSmNT,nSample)),
                                  legacy = TRUE)
            }

            GDCdownload(queryDown, directory = DataDirectory)

            RSEobject <- GDCprepare(query = queryDown, 
                                   save = TRUE, 
                                   directory =  DataDirectory,
                                   save.filename = FileNameData)
        }else{
            load(FileNameData)
            data <- get("data")
            RSEobject <- data
        }

        dataPrep <- TCGAanalyze_Preprocessing(object = RSEobject,
                                              cor.cut = cor.cut,
                                              datatype = "raw_count")                      

        dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                              geneInfo = geneInfo,
                                              method = "gcContent")                

        dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                          method = "quantile", 
                                          qnt.cut =  qnt.cut)  

        return(dataFilt)

    }else if(dataType == "Methylation"){
        FileNameData <- paste0(DataDirectory, "_","Illumina Human Methylation 27",".rda")
        if(!file.exists(FileNameData)){
            query.met <- GDCquery(project = CancerProject, 
                                legacy = TRUE,
                                data.category = "DNA methylation",
                                platform = "Illumina Human Methylation 27")

            query.met.samples <-  query.met$results[[1]]$cases

            dataSmTP <- TCGAquery_SampleTypes(query.met.samples,"TP")
            dataSmNT <- TCGAquery_SampleTypes(query.met.samples,"NT")

            message <- paste0(CancerProject," ... with ", length(dataSmTP), " TP, ", length(dataSmNT), " NT")
            print(message)


            if(is.null(nSample)){
                sampDown <- c(dataSmTP, dataSmNT)
            }else{
                sampDown <- c(sample(dataSmTP,nSample), sample(dataSmNT,nSample))
            }

            query.met <- GDCquery(project = CancerProject, 
                                legacy = TRUE,
                                data.category = "DNA methylation",
                                platform = "Illumina Human Methylation 27",
                                barcode = sampDown)

            GDCdownload(query.met)

            cancer.met <- GDCprepare(query = query.met,
                                 save = TRUE, 
                                 save.filename = FileNameData,
                                 summarizedExperiment = TRUE)
                                 
            return(cancer.met)
      
        }else{
            load(FileNameData)
            RSEobject <- data
        }
    }
}