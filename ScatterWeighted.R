
#4)fem el subsets per families i per adj Pval
scatter3.PL <- data.scatter3[which(data.scatter3$FamShort == "PL"),]


#5)Fem els gràfics per separat
p.PL <- ggplot(scatter3.PL, aes(x = NumberOfCarbons, y = NumberOfDoubleBonds, 
                              size = abs(logFC), fill = adj.P.Val)) +
      geom_point(shape = 21) +
      ggtitle("Phospholipids") +
      labs(x = "No. of Carbons") +
      labs(size = "log Fold Change", fill = "Adjusted pValue") +
      scale_fill_continuous(high = "#f2e2ff", low = "#ab56ef") + 
      scale_size_continuous(range = c(1,5)) +
      scale_y_continuous("No. of Double Bonds",limits = c(0,7)) +
      theme_bw() 

pdf(file.path(resultsDir,"WeightedScatter.pdf"))
grid.arrange(p.PL, p.ST, p.DAG, p.TAG, p.SM, ncol = 2)
dev.off()
