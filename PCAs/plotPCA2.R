#' Plot PCA
#'
#' This function plot the two first principal components of a dataset
#'
#' @param datos A data set where only the first column is character (usually sample names)
#' @param labels A vector with the sample names (usually in targets$shortname column)
#' @param factor A vector with the different levels to study in the PCA (usually in targets$Grupo column)
#' @param title The name of the factor which you are colouring the sample
#' @param scale The data needs to be scaled? Values= TRUE/FALSE
#' @export plotPCA2
#' @import ggplot2 ggrepel
#' @examples
#' plotPCA2(datos=ExampleData$data,labels = ExampleData$targets$shortname,
#' factor=ExampleData$targets$grupo,title = "Grupo",scale=FALSE)

plotPCA2 <- function ( datos, labels, factor,title,scale) {
  data <- prcomp(t(datos[,-1]),scale=scale)
  #ajustos del gràfic
  dataDf <- data.frame(data$x)
  Grupo <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  # the graphic
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Grupo), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5))
  # the graphic with ggrepel
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels)) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("Principal Component Analysis for:",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5))
}
