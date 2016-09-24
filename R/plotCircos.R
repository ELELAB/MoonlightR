#' @title plotCircos
#' @description
#' This function visualize the plotCircos
#' @param listMoonlight output Moonlight function 
#' @param additionalFilename additionalFilename 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize circos.par
#' @importFrom circlize circlize
#' @importFrom circlize circos.clear
#' @importFrom circlize circos.initialize
#' @importFrom circlize circos.trackPlotRegion get.cell.meta.data circos.text circos.rect
#' @importFrom grDevices dev.off pdf rainbow rgb
#' @export
#' @return no return value, plot is saved
#' @examples 
#' plotCircos(listMoonlight = listMoonlight, additionalFilename = "_ncancer5")

plotCircos <- function(listMoonlight, additionalFilename = NULL){

    ### prepare data
    n <- length(listMoonlight)
    mycancertypes <- names(listMoonlight)

    # listMoonlight$listMoonlight.tsg <- listMoonlight[[1]]
    # listMoonlight$listMoonlight.osg <- listMoonlight[[2]]

    mytsg <- myocg <- NULL
    for(i in 1:n){
        mytsg <- c(mytsg, list(names(listMoonlight[[i]]$listCandidates$TSG)))
        myocg <- c(myocg, list(names(listMoonlight[[i]]$listCandidates$OCG)))

    }
    # mytsg <-  listMoonlight[[2]]
    # myocg <-  listMoonlight[[1]]

    n.mygenes <- sapply(mytsg, length) + sapply(myocg,length)
    mynames <- mycancertypes
    ind.rm <- which(n.mygenes==0)

    if(length(ind.rm)>0){
        mytsg <-  listMoonlight[[2]][-ind.rm]
        myocg <-  listMoonlight[[1]][-ind.rm]
        n.mygenes <- sapply(mytsg, length) + sapply(myocg,length)

        mynames <- mycancertypes[-ind.rm]
    }
    n <- n - length(ind.rm)

    ntsg <- sapply(mytsg,length)
    nocg <- sapply(myocg, length)


    pdf(paste0("circos_ocg_tsg",additionalFilename,".pdf"), width=16, height=16)
    df1 <- data.frame("order" =c(1:n), "cancertype"=mynames, "xmin"=rep(0,n), "xmax"=n.mygenes)

    ### Plot sectors (outer part)
    par(mar=rep(0,4))
    circlize::circos.clear()

    ### Basic circos graphic parameters
    circlize::circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)

    if(n>12){
        mycols <- rainbow(n+4)[1:n]
    }else{
        mycols <- RColorBrewer::brewer.pal(n, "Set3")
    }
    ### Sector details
    circlize::circos.initialize(factors = df1$cancertype, xlim = cbind(df1$xmin, df1$xmax))

    ### Plot sectors
    circlize::circos.trackPlotRegion(ylim = c(0, 1), factors = df1$cancertype, track.height=0.1,
                      #panel.fun for each sector
                      panel.fun = function(x, y) {
                      #select details of current sector
                      name = get.cell.meta.data("sector.index")
                      i = get.cell.meta.data("sector.numeric.index")
                      xlim = get.cell.meta.data("xlim")
                      ylim = get.cell.meta.data("ylim")

                      #text direction (dd) and adjusmtents (aa)
                      theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                      dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                      aa = c(1, 0.5)
                      if(theta < 90 || theta > 270)  aa = c(0, 0.5)

                      circlize::circos.text(x=mean(xlim), y=1.7, labels=paste0(name,"(",sapply(myocg,length)[i],", ",sapply(mytsg,length)[i],")"), facing = dd, cex=0.75,  adj = aa)

                      #plot main sector
                      print(df1$rcol[i])
                      circlize::circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2],
                                  col = mycols[i], border=mycols[i])


                      circlize::circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-sapply(mytsg, length)[i], ytop=ylim[1]+0.3, 
                                  col = "darkgreen", border = "darkgreen")
                      circlize::circos.rect(xleft=sapply(myocg, length)[i], ybottom=ylim[1], xright=xlim[2], ytop=ylim[1]+0.3, 
                                  col = "goldenrod", border = "goldenrod")

                      # #white line all the way around
                      circlize::circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")

                      
                    })

    ## plot inner part
    ### OCG - OCG
    mycol.ocg <- rgb(10/255, 99/255, 12/255, 0.5)
    for(i in 1:(n-1)){
        print(paste("cancertype", i, "out of", n))
        for(j in (i+1):n){
            if(i!=j){
                ind <- which(myocg[[i]] %in% myocg[[j]])
                for(k in ind){
                    ind2 <- which(myocg[[j]]==myocg[[i]][k])
                    circlize::circos.link(sector.index1=df1$cancertype[i], point1=c(k-1,k),
                        sector.index2=df1$cancertype[j], point2=c(ind2-1,ind2), col = mycol.ocg)
                }
            }
        }
    }
    
    ### TSG - TSG 
    mycol.tsg <- rgb(217/255, 164/255, 50/255, 0.5)
    for(i in 1:(n-1)){
        print(paste("cancertype", i, "out of", n))
        for(j in (i+1):n){
            if(i!=j){
                ind <- which(mytsg[[i]] %in% mytsg[[j]])
                for(k in ind){
                    ind2 <- which(mytsg[[j]]==mytsg[[i]][k])
                    circlize::circos.link(sector.index1=df1$cancertype[i], point1=c(sapply(myocg, length)[i]+k-1, sapply(myocg, length)[i]+ k),
                        sector.index2=df1$cancertype[j], point2=c(sapply(myocg, length)[j]+ind2-1,sapply(myocg, length)[j]+ind2), col = mycol.tsg)
                }
            }
        }
    }

    # ### TSG - OSG
    mycol.tsg.osg <- rgb(252/255, 51/255, 57/255, 0.5)
    for(i in 1:n){
        print(paste("cancertype", i, "out of", n))
        for(j in 1:n){
            ind <- which(myocg[[i]] %in% mytsg[[j]])
            for(k in ind){
                ind2 <- which(mytsg[[j]]==myocg[[i]][k])
                circlize::circos.link(sector.index1=df1$cancertype[i], point1=c(k-1,k),
                    sector.index2=df1$cancertype[j], point2=c(sapply(myocg, length)[j]+ind2-1,sapply(myocg, length)[j]+ind2), col = mycol.tsg.osg)
            }
        }
    }

    dev.off()
}
