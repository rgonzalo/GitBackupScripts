
#datos = no se resta la primera columna
#labels = nom de les mostres
#factor = factor por el que se colorearán las muestras
#tittle = pequeño título que se añade a "Principal Component Analysis for: "
#glineas = se controla el grosor de las líneas de las etiquetas (segment.size)
#scale = si quieres escalar los datos o no
#size = se controla el tamaño de las etiquetas (size=1.5)
#colores <- vector con el color de cada factor de forma ordenada
require(ggplot2)
require(ggrepel)

plotPCA3 <- function (datos, labels, factor,title,scale,colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos),scale=scale)
  #ajustos del gràfic
  dataDf <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  # the graphic
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5))
  # the graphic with ggrepel
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25,size=1.5) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("Principal Component Analysis for: ",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colores)
  }
