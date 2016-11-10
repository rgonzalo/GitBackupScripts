library(heatmaply)

#load the normalized data (example from TGarcia miRNAs study)
norm_data<-read.table("FF_rma_anotations.csv",sep="\t",header=TRUE)
colnames(norm_data)
rownames(norm_data)<-norm_data$Probe.Set.ID ; norm_data<-norm_data[,-c(1,14:27)]

#load the targets file
my.targets <-read.table("targets1.csv",header = TRUE, row.names = 1, sep="\t") 
colnames(norm_data)<-rownames(my.targets)


#fit contrast to select genes by fold change or significance
require(Biobase)
stopifnot(require(limma))
factArea <- as.factor(my.targets$Area)
factID <-as.factor(my.targets$Patient)
design <- model.matrix(~0+factID+factArea)

rownames(design)<-colnames(norm_data)
print(design)

#CONTRASTES (no se ha creado objeto eset y se han puesto datos normalizados a saco)
require(limma)
comparison<-"FF"
cont.matrix1 <- makeContrasts (ICvsCL=factAreaIC,
                               levels=design)
fit1<-lmFit(norm_data, design)
fit.main1<-contrasts.fit(fit1, cont.matrix1)
fit.main1<-eBayes(fit.main1)

topTab <-  topTable (fit.main1, number=nrow(fit.main1), adjust="fdr")

#select the genes to depict in the heatmap
selected<-rownames(topTab[topTab$P.Value<0.05 & topTab$logFC>abs(0.5),])
length(selected)#32
exprs2cluster1 <-as.matrix(norm_data)[selected,]

#draw the heatmap
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)
heatmaply(exprs2cluster1,my_palette)
