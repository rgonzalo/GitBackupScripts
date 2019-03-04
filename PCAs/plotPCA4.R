require(ggplot2)

plotPCA4 <- function (datos, factor, title, scale, colores) {
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
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  p1 + labs(x = paste0("PC1 ",loads[1],"%"), y = paste0("PC2 ",loads[2],"%")) +
    scale_color_manual(values=colores) +
    ggtitle(paste("Principal Component Analysis for: ",title,sep=" ")) +
    theme(plot.title = element_text(hjust = 0.5)) 
}

