require(ggplot2)

#PCA without labels
plotPCA4 <- function ( datos, labels, factor,title,scale) {
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
  p1 + labs(x = paste0("PC1 ",loads[1],"%"), y = paste0("PC2 ",loads[2],"%"))   }