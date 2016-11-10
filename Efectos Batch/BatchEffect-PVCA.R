################################
########     PVCA       ########
################################
### Determinar fonts de variació,
### quantificant la proporció de 
### variabilitat explicada per
### cada variable
################################
## http://www.bioconductor.org/packages/release/bioc/manuals/pvca/man/pvca.pdf
################################

if (!require(pvca)) install.packages("pvca") 
load("dadesexemple/normalizedData2014.Rda")
targets2014 <- read.csv2 ( "dadesexemple/targets2014.csv",sep="\t")
pData(my.norm) <- targets2014 
pct_threshold <- 0.6
batch.factors <- c("Grupo", "Batch", "Gender","Ficoll")
pvcaObj <- pvcaBatchAssess (my.norm, batch.factors, pct_threshold)

bp <- barplot(pvcaObj$dat, xlab = "Effects",
              ylab = "Weighted average proportion variance",
              ylim= c(0,1.1),col = c("blue"), las=2,
              main="PVCA estimation bar chart (2014)")
axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.7, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)
