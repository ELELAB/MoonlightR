#' @title Pattern Recognition Analysis (PRA)
#' @description
#' This function carries out the pattern recognition analysis
#' @param dataURA output URA function
#' @param BPname BPname
#' @param thres.role thres.role
#' @param seed seed value
#' @importFrom randomForest randomForest
#' @export
#' @return returns list of TSGs and OCGs when biological processes are provided, otherwise a randomForest based classifier that can be used on new data
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
PRA <- function(dataURA, BPname, ntree = 1000, thres.role = 0, seed=12345){
    set.seed(seed)

    names.blocking <- tabGrowBlock[which(tabGrowBlock$Cancer.blocking == "Increasing"), "Disease"]
    names.growing <- tabGrowBlock[which(tabGrowBlock$Cancer.growing == "Increasing"), "Disease"]

    if(is.null(BPname)){

      	# print("random forest")
        common.genes.tsg <- intersect(knownDriverGenes$TSG, rownames(dataURA))
        common.genes.ocg <- intersect(knownDriverGenes$OCG, rownames(dataURA))

        dataTrain <- data.frame(dataURA[c(common.genes.tsg, common.genes.ocg),],
            "labels"=as.factor(c(rep(1,length(common.genes.tsg)),rep(0,length(common.genes.ocg)))))

        fit.rf <- randomForest((labels)~.,data=dataTrain, importance=TRUE, ntree = ntree)
        return(fit.rf)

    }else{
    	# print("selected processes")
        ind.proc.growing <- which(BPname %in% names.growing)
        ind.proc.blocking <- which(BPname %in% names.blocking)

        names.genes.tsg <- names(which(dataURA[,BPname[ind.proc.growing]] < -thres.role & dataURA[,BPname[ind.proc.blocking]] > thres.role))
        names.genes.ocg <- names(which(dataURA[,BPname[ind.proc.growing]] > thres.role & dataURA[,BPname[ind.proc.blocking]] < -thres.role))

        if(length(names.genes.tsg)>0){
            scores.tsg <- apply(abs(dataURA[names.genes.tsg,]), 1, mean)
        }else{
            scores.tsg <- NULL
        }
        
        if(length(names.genes.ocg)>0){
            scores.ocg <- apply(abs(dataURA[names.genes.ocg,]), 1, mean)   
        }else{
            scores.ocg <- 0
        }
        
        return(list("TSG"=scores.tsg, "OCG"=scores.ocg))
    	# dataURA[,BPname]
    	
    }


}