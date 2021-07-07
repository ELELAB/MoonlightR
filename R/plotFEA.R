#' @title plotFEA
#' @description
#' This function visualize the functional enrichment analysis (FEA)'s barplot
#' @param dataFEA dataFEA
#' @param topBP topBP
#' @param additionalFilename additionalFilename
#' @param height Figure height
#' @param width Figure width
#' @param offsetValue offsetValue
#' @param angle angle
#' @param xleg xleg
#' @param yleg yleg
#' @param titleMain title of the plot
#' @param minY minY
#' @param maxY maxY
#' @param mycols colors to use for the plot
#' @importFrom graphics barplot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics text
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export
#' @return no return value, FEA result is plotted
#' @examples
#' dataFEA <- FEA(DEGsmatrix = DEGsmatrix)
#' plotFEA(dataFEA = dataFEA, additionalFilename = "_example",height = 20,width = 10)
plotFEA <- function (dataFEA, topBP = 10, additionalFilename = NULL, height,
                     width, offsetValue = 5, angle = 90, xleg = 35, yleg = 5, titleMain,
                     minY = -5, maxY = 10,  mycols = c("#8DD3C7", "#FFFFB3", "#BEBADA")) {
    titleMain <- "TCGA BRCA DEGs"
    if (!is.null(additionalFilename)) {
        pdf(additionalFilename, width, height)
    }
    par(mar = c(12, 5, 5, 1))
    dataFEA <- dataFEA[order(abs(dataFEA$Moonlight.Z.score),
                             decreasing = TRUE), ]
    tmp <- dataFEA[1:topBP, ]
    tmp <- as.data.frame(tmp)
    tmp$FDR <- as.numeric(tmp$FDR)
    tmp$Moonlight.Z.score <- as.numeric(tmp$Moonlight.Z.score)
    tmp$commonNg <- as.numeric(tmp$commonNg)
    tmp$FDR <- -log2(tmp$FDR)/10
    tmp <- tmp[order(tmp[,4], decreasing = TRUE),]
    toPlot <- matrix(0, nrow = 3, ncol = nrow(tmp))
    toPlot[1, ] <- tmp$FDR
    toPlot <- toPlot[, order(toPlot[1, ], decreasing = TRUE)]
    toPlot[2, as.numeric(tmp$Moonlight.Z.score) > 0] <- as.numeric(tmp$Moonlight.Z.score[as.numeric(tmp$Moonlight.Z.score) >
                                                                                             0])
    toPlot[3, as.numeric(tmp$Moonlight.Z.score) < 0] <- as.numeric(tmp$Moonlight.Z.score[as.numeric(tmp$Moonlight.Z.score) <
                                                                                             0])
    toPlot[1, which(toPlot[1, ] > 10)] <- 10
    toPlot <- as.matrix(toPlot)
    xAxis <- barplot(toPlot, beside = TRUE, col = mycols, ylab = "-log2FDR/10 and Moonlight z-score",
                     names = NULL, main = paste0("FEA - Enriched BioFunctions", titleMain,
                                                 sep = " "), ylim = c(minY, maxY))
    legend("topright", bty = "n", xleg, yleg, legend = c("-log2FDR/10",
                                                         "Moonlight Z-score Increased", "Moonlight Z-score Decreased"),
           fill = mycols, cex = 0.5)
    text(xAxis[2, ], par("usr")[3] + offsetValue, srt = angle,
         adj = 1, xpd = TRUE, labels = paste(tmp$Diseases.or.Functions.Annotation,
                                             " (n=", tmp$commonNg, ")", sep = ""), cex = 1)
    if (!is.null(additionalFilename)) {
        dev.off()
    }
}
