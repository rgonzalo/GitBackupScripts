library(VennDiagram)#per fer un Venn una mica més bonic
library(limma)#venn clàsic


#Exemple per fer-ho amb tres listes a partir csv creats amb el Basic Pipe
################################################################################
top1 <- read.csv2(file = "NOM DE L'ARXIU1.csv",header = TRUE, sep = ";", dec=",")
top2 <- read.csv2(file = "NOM DE L'ARXIU2.csv",header = TRUE, sep = ";", dec=",")
top3 <- read.csv2(file = "NOM DE L'ARXIU3.csv",header = TRUE, sep = ";", dec=",")

dim(top1)
dim(top2)
dim(top3)
head(top1)
head(top2)
head(top3)

####seleccionem els cutoffs
######################################################################
pass1 <- top1[which(top1$P.Value<0.01),]
tail(pass1)
pass1 <- pass1[which(pass1$logFC > 1 | pass1$logFC < -1),]
tail(pass1)

pass2 <- top2[which(top2$P.Value<0.01),]
tail(pass2)
pass2 <- pass2[which(pass2$logFC > 1 | pass2$logFC < -1),]
tail(pass2)

pass3 <- top3[which(top3$P.Value<0.01),]
tail(pass3)
pass3 <- pass3[which(pass3$logFC > 1 | pass3$logFC < -1),]
tail(pass3)

####mirem que no entrades repetides per Affy ID i s'ajunten
######################################################################
list1 <- as.character(sort(unique(pass1$X)))
length(list1)
list2 <- as.character(sort(unique(pass2$X)))
length(list2)
list3 <- as.character(sort(unique(pass3$X)))
length(list3)

list <- c(list1, list2, list3)
length(list)
list <- sort(unique(list))
length(list)

####es crea un data.frame que omplim de 0 i després de 1 si hi coexisteixen en les dues llistes
###############################################################################################

df <- data.frame(micros = list, l1 = rep(0,length(list)), l2 = rep(0,length(list)), l3 = rep(0,length(list)))
head(df)

df$l1 <- as.numeric((df$micros %in% list1)==T)
df$l2 <- as.numeric((df$micros %in% list2)==T)
df$l3 <- as.numeric((df$micros %in% list3)==T)

####es fa el Venn i s'exporta a un arxiu
#############################################################################
vennDiagram(vennCounts(df[,2:4]), names=c("Comp1","Comp2","Comp3"),
            circle.col=c("red", "green", "yellow"))

pdf("VennDiagram.pdf")
vennDiagram(vennCounts(df[,2:4]), names=c("Comp1","Comp2","Comp3"),
            circle.col=c("red", "green", "yellow"))
dev.off()

##així es faria utilitzant la funcio del paquet VennDiagram (OJO ESTA ES PARA DOS LISTAS)
#############################################################################
overlap<-calculate.overlap(x=list("list1"=list1,"list2"=list2))

draw.pairwise.venn(length(overlap$a1),length(overlap$a2),length(overlap$a3),
                   category=c("Comp1","Comp2"),scaled = TRUE,euler.d = TRUE, 
                   fill = c("blue", "red"),lty = "blank",cat.pos = c(190, 0))

pdf("VennDiagram.pdf")
draw.pairwise.venn(length(overlap$a1),length(overlap$a2),length(overlap$a3),
                   category=c("Comp1","Comp2"),scaled = TRUE,euler.d = TRUE, 
                   fill = c("blue", "red"),lty = "blank",cat.pos = c(190, 0))
dev.off()

##es grava l'arxiu de text on s'indica quin gen hi és a cada venn
#############################################################################
datos<-top2[which(top2$X %in% df$AffyID),]
datos<-datos[,-c(3:length(colnames(datos)))] 
colnames(datos)<-c("AffyID","Symbols") 
datos<-merge(datos,df,by="AffyID") 
colnames(datos)<-c("AffyID","Symbols","Comp1","Comp2")
write.csv(datos, file="venn_genes.csv")
             