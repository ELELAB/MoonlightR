#' MoonlightR
#'
#' @name MoonlightR
#' @docType package
NULL


#' Information on 101 biological processes
#'
#' A data set containing the following data: 
#'
#' \itemize{
#'		\item DiseaseList list for 101 biological processes, each containing a matrix with five columns: ID, Genes.in.dataset, Prediction based on expression direction, Log ratio, Findings
#' }
#' @docType data
#' @keywords datasets
#' @name DiseaseList
#' @usage data(DiseaseList)
#' @format A list of 101 matrices
#' @return list of 101 matrices
NULL

#' DEG Differentially expressed genes
#'
#' A data set containing the following data: 
#'
#' \itemize{
#'		\item DEGsmatrix matrix with 3502 rows (genes) and five columns  "logFC"  "logCPM" "LR"     "PValue" "FDR" 
#' }
#' @docType data
#' @keywords datasets
#' @name DEGsmatrix
#' @usage data(DEGsmatrix)
#' @format A 3502x5 matrix
#' @return the 3502x5 matrix
NULL

#' Gene Expression (Rnaseqv2) data from TCGA LUAD 
#'
#' A data set containing the following data: 
#'
#' \itemize{
#'		\item dataFilt matrix with 13742 rows (genes) and 20 columns samples with TCGA's barcodes (10TP, 10NT)
#' }
#' @docType data
#' @keywords datasets
#' @name dataFilt
#' @usage data(dataFilt)
#' @format A 13742x20 matrix
#' @return a 13742x20 matrix
NULL

#' Output list from Moonlight
#'
#' A list containing the following data: 
#'
#' \itemize{
#'		\item listMoonlight output from moonlight's pipeline containing dataDEGs, dataURA, listCandidates
#' }
#' @docType data
#' @keywords datasets
#' @name listMoonlight
#' @usage data(listMoonlight)
#' @format A Large list with 5 elements
#' @return output from moonlight pipeline
NULL


#' Information growing/blocking characteristics for 101 selected biological processes
#'
#' A data set containing the following data: 
#'
#' \itemize{
#'	\item tabGrowBlock matrix that defines if a process is growing or blocking cancer development, for each 101 biological processing 
#' }
#' @docType data
#' @keywords datasets
#' @name tabGrowBlock
#' @usage data(tabGrowBlock)
#' @format A 101x3 matrix
#' @return a 101x3 matrix
NULL

#' Information on GEO data (and overlap with TCGA)#'
#' A data set containing the following data: 
#'
#' \itemize{
#'	\item GEO_TCGAtab a 18x12 matrix that provides the GEO data set we matched to one of the 18 given TCGA cancer types 
#' }
#' @docType data
#' @keywords datasets
#' @name GEO_TCGAtab
#' @usage data(GEO_TCGAtab)
#' @format A 101x3 matrix
#' @return a 101x3 matrix
NULL

#' GRN gene regulatory network output
#'
#' output from GRN function 
#'
#' \itemize{
#'		\item dataGRN list of 2 elements miTFGenes, maxmi from GRN function
#' }
#' @docType data
#' @keywords datasets
#' @name dataGRN
#' @usage data(dataGRN)
#' @format A large list of 2 elements 
#' @return a large list of 2 elements
NULL


#' Information on known cancer driver gene from COSMIC
#'
#' A data set containing the following data: 
#'
#' \itemize{
#'	\item TSG known tumor suppressor genes
#'  \item OCG known oncogenes
#' }
#' @docType data
#' @keywords datasets
#' @name knownDriverGenes
#' @usage data(knownDriverGenes)
#' @format A 101x3 matrix
#' @return a 101x3 matrix
NULL


#' Information about genes
#'
#' A data set containing the following data: 
#'
#' \itemize{
#'		\item EAGenes matrix with 20038 rows (genes) and five columns   
#'		"ID"          "Gene"        "Description" "Location"    "Family"     
#' }
#' @docType data
#' @keywords datasets
#' @name EAGenes
#' @usage data(EAGenes)
#' @format A 20038x5 matrix
#' @return a 20038x5 matrix
NULL

#' Information about genes for normalization
#'
#' A data set containing the following data: 
#'
#' \itemize{
#'		\item geneInfo matrix with 20531 rows (genes) and 3 columns   
#'		"geneLength" "gcContent"  "chr"         
#' }
#' @docType data
#' @keywords datasets
#' @name geneInfo
#' @usage data(geneInfo)
#' @format A data frame with 20531 rows and 3 variables
#' @return a 20531x3 matrix
NULL

#' Information on GDC projects
#'
#' A character vector of GDC projects: 
#'
#' \itemize{
#'		\item character vector for GDC projects.
#' }
#' @docType data
#' @keywords datasets
#' @name GDCprojects
#' @usage data(GDCprojects)
#' @format A character vector of 39 elements
#' @return character vector of 39 elements
NULL

#' Output example from function Upstram Regulator Analysis
#'
#' A data set containing the following data: 
#'
#' \itemize{
#'		\item dataURA matrix with 100 rows (genes) and 2 columns   
#'		 "apoptosis"              "proliferation of cells"
#' }
#' @docType data
#' @keywords datasets
#' @name dataURA
#' @usage data(dataURA)
#' @format A data frame with 100 rows and 2 variables
#' @return a 100x2 matrix
NULL