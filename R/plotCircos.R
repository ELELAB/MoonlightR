#' @title plotCircos
#' @description
#' This function visualize the plotCircos
#' @param listMoonlight output Moonlight function 
#' @param additionalFilename additionalFilename 
#' @param intensityColOCG intensityColOCG
#' @param intensityColTSG intensityColTSG
#' @param intensityColDual intensityColDual
#' @param fontSize fontSize
#' @param listMutation listMutation    
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize circos.par
#' @importFrom circlize circlize
#' @importFrom circlize circos.clear
#' @importFrom circlize circos.initialize
#' @importFrom circlize circos.trackPlotRegion get.cell.meta.data circos.text circos.rect
#' @importFrom grDevices dev.off rainbow rgb
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export
#' @return no return value, plot is saved
#' @examples 
#' plotCircos(listMoonlight = listMoonlight, additionalFilename = "_ncancer5")

plotCircos <- function(listMoonlight, listMutation = NULL, additionalFilename = NULL, 
                       intensityColOCG = 0.5, intensityColTSG = 0.5, intensityColDual = 0.5, 
                       fontSize=1){

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
    names(mytsg) <- names(myocg) <- mycancertypes
    # mytsg <-  listMoonlight[[2]]
    # myocg <-  listMoonlight[[1]]

    n.mygenes <- sapply(mytsg, length) + sapply(myocg,length)
    mynames <- mycancertypes
    ind.rm <- which(n.mygenes==0)

    if(length(ind.rm)>0){
        mytsg <-  mytsg[-ind.rm]
        myocg <-  myocg[-ind.rm]
        n.mygenes <- sapply(mytsg, length) + sapply(myocg,length)

        mynames <- mycancertypes[-ind.rm]
        listMutation <- listMutation[-ind.rm]
    }
    n <- n - length(ind.rm)

    ntsg <- sapply(mytsg,length)
    nocg <- sapply(myocg, length)


    pdf(paste0("circos_ocg_tsg",additionalFilename,".pdf"), width=16, height=16)
    df1 <- data.frame("order" =c(1:n), "cancertype"=mynames, "xmin"=rep(0,n), "xmax"=n.mygenes)

    ### Plot sectors (outer part)
    par(mar=c(1,1,6,6))
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

                        circlize::circos.text(x=mean(xlim), y=1.7, labels=paste0(name,"\n(",sapply(myocg,length)[i],", ",sapply(mytsg,length)[i],")"), facing = dd, cex=fontSize,  adj = aa)

                        #plot main sector
                        # print(df1$rcol[i])
                        print(xlim)
                        circlize::circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2],
                                    col = mycols[i], border=mycols[i])


                        circlize::circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-sapply(mytsg, length)[i], ytop=ylim[1]+0.3, 
                                    col = "darkgreen", border = "darkgreen")
                        circlize::circos.rect(xleft=sapply(myocg, length)[i], ybottom=ylim[1], xright=xlim[2], ytop=ylim[1]+0.3, 
                                    col = "goldenrod", border = "goldenrod")

                        # #white line all the way around
                        circlize::circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")



                      
                    })  
    # for(i in 1:n){
    if(!is.null(listMutation)){
      mycols.mut <- c("white","darkviolet","tomato3")
      for(i in 1:length(listMutation)){
        print(i)
        print(dim(listMutation[[i]]))
        if(length(myocg[[i]])>0){
            myoverlap <-  myocg[[i]] %in% unique(as.matrix(listMutation[[i]][which(listMutation[[i]][,"Consequence"]=="inframe_deletion"),1]))
            myoverlap2 <-  myocg[[i]] %in% unique(as.matrix(listMutation[[i]][which(listMutation[[i]][,"Consequence"]=="inframe_insertion"),1]))
            myoverlap3 <-  myocg[[i]] %in% unique(as.matrix(listMutation[[i]][which(listMutation[[i]][,"Consequence"]=="missense_variant"),1]))

            cnt <- 0
            print(length(myocg[[i]]))
            print(length(myoverlap))
            print(length(myoverlap2))

            for(j in 1:length(myocg[[i]])){
                # print(j)
                if(myoverlap[j]){
                  circlize::circos.rect(sector.index = df1$cancertype[i], xleft=cnt, ybottom=-1/2, xright=cnt+1, ytop=-0.33, 
                                        col = mycols.mut[2], border = mycols.mut[2])
                }
                if(myoverlap2[j]){
              circlize::circos.rect(sector.index = df1$cancertype[i],xleft=cnt, ybottom=-0.77, xright=cnt+1, ytop=-0.6, 
                                        col = mycols.mut[2], border = mycols.mut[2])
            }
           if(myoverlap3[j]){
              circlize::circos.rect(sector.index = df1$cancertype[i],xleft=cnt, ybottom=-1.04, xright=cnt+1, ytop=-0.87, 
                                        col = mycols.mut[2], border = mycols.mut[2])
          }

              cnt <- cnt+1
            }
        }

        if(length(mytsg[[i]])>0){
            myoverlap <-  mytsg[[i]] %in% unique(as.matrix(listMutation[[i]][which(listMutation[[i]][,"Consequence"]=="inframe_deletion"),1]))
            myoverlap2 <-  mytsg[[i]] %in% unique(as.matrix(listMutation[[i]][which(listMutation[[i]][,"Consequence"]=="inframe_insertion"),1]))
            myoverlap3 <-  mytsg[[i]] %in% unique(as.matrix(listMutation[[i]][which(listMutation[[i]][,"Consequence"]=="missense_variant"),1]))

            cnt <- 0
            print(length(mytsg[[i]]))
            print(length(myoverlap))
            print(length(myoverlap2))

            for(j in 1:length(mytsg[[i]])){
                # print(j)
                if(myoverlap[j]){
                  circlize::circos.rect(sector.index = df1$cancertype[i], xleft=nocg[i]+cnt, ybottom=-1/2, xright=nocg[i]+cnt+1, ytop=-0.33, 
                                        col = mycols.mut[3], border = mycols.mut[3])
                }
                if(myoverlap2[j]){
              circlize::circos.rect(sector.index = df1$cancertype[i],xleft=nocg[i]+cnt, ybottom=-0.77, xright=nocg[i]+cnt+1, ytop=-0.6, 
                                        col = mycols.mut[3], border = mycols.mut[3])
            }
           if(myoverlap3[j]){
              circlize::circos.rect(sector.index = df1$cancertype[i],xleft=nocg[i]+cnt, ybottom=-1.04, xright=nocg[i]+cnt+1, ytop=-0.87, 
                                        col = mycols.mut[3], border = mycols.mut[3])
          }

              cnt <- cnt+1
            }
        }

        
    }
    }

    ## plot inner part
    ### OCG - OCG
    mycol.ocg <- rgb(10/255, 99/255, 12/255, intensityColOCG)
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
    mycol.tsg <- rgb(217/255, 164/255, 50/255, intensityColTSG)
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
    mycol.tsg.osg <- rgb(252/255, 51/255, 57/255, intensityColDual)
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

    # if( (which = dev.cur()) != 1){graphics.off() }
    dev.off()

    }
