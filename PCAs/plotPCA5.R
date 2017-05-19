#PCA donde se puede seleccionar los colores de los puntos:scale_color_manual
#colores <- vector con el color de cada factor de forma ordenada
require(ggplot2)
require(ggrepel)

plotPCA5 <- function (datos, labels, factor,title,scale,colores) {
  data <- prcomp(t(datos),scale=scale)
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
    coord_cartesian(xlim = c(min(data$x[,1])-1,max(data$x[,1])+1))
  # the graphic with ggrepel
  p1 + geom_text_repel(aes(y = PC2 + 0, label = labels),segment.size = 0.25,size=2.5) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("Principal Component Analysis for:",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colores)
}