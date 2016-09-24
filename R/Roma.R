#' Roma
#'
#' This function carries out the ROMA analysis
#' @param exp_data expression data where the activity of the signatures needs to be evaluated
#' @param signature in gmt format
#' @param path absolute path to MoonlightR folder
#' @export 
#' @return matrix of roma scores
#' @examples 
#' dataDEGs <- DPA(dataFilt = dataFilt, dataType = "Gene expression")
#' # to change with ROMA

Roma <- function(exp_data,signature,path){
    Temp_fold<-tempdir()
    write.table(exp_data,paste(Temp_fold,'exp.txt',sep=''),sep="    ",row.names=TRUE,col.names=TRUE,quote=FALSE)
    write.table(signature,paste(Temp_fold,'signature.gmt',sep=''),sep="    ",row.names=TRUE,col.names=TRUE,quote=FALSE)
    system(paste("export CLASSPATH=",path,"Roma-master/jar/ROMA.jar:",path,"Roma-master/lib/VDAOEngine.jar",sep=""))
    system(paste("java -Xmx8024m -cp ",path,"Roma-master/jar/ROMA.jar fr.curie.ROMA.ModuleActivityAnalysis -dataFile ",Temp_fold,"exp.txt -moduleFile ",Temp_fold,"signature.gmt -outputFolder ",Temp_fold,"out/ -centerData 1 -numberOfPermutations 1000",sep=""), intern=TRUE)
    out<-as.matrix(read.table(paste(Temp_fold,"out/moduletable_withscores.txt",sep=""),header=TRUE,sep="    ",row.names=1))
    L1<-as.matrix(out$L1)
    rownames(L1)<-rownames(out)
    return(L1)
}