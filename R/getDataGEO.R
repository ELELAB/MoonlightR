#' getDataGEO
#'
#' This function retrieves and prepares GEO data
#' @param GEOobject GEOobject 
#' @param platform platform
#' @param TCGAtumor tumor name
#' @importFrom GEOquery getGEO
#' @export
#' @return return GEO gset
#' @examples
#' dataGEO <-  getDataGEO(GEOobject = "GSE20347",platform = "GPL571")

getDataGEO <- function(GEOobject = "GSE39004", platform = "GPL6244", TCGAtumor=NULL){
    
  GEO_TCGAtab <- get("GEO_TCGAtab")
  
    if (length(TCGAtumor)!=0){
        GEOobject <- GEO_TCGAtab[GEO_TCGAtab$Cancer ==  TCGAtumor,"Dataset"]
        platform <- GEO_TCGAtab[GEO_TCGAtab$Cancer ==  TCGAtumor,"Platform"]
    }

    gset <- getGEO(GEOobject, GSEMatrix =TRUE, AnnotGPL=TRUE)
    if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    return(gset) 
}