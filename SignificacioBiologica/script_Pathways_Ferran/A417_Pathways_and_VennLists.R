###################################################
###################################################
##step 0: Setup
###################################################

## set working directory wherever you like
#setwd("/media/ferran/Datos/estudis/RNAseq/2015-09-XavierSantamaria-IVI-A319")

workingDir <- getwd() # Compte, cal fixar-lo abans
resultsDir <- file.path(workingDir, "results")

#source("http://bioconductor.org/biocLite.R")



###################################################
###################################################
##step 1: READ EXPRESSIONS FROM csv FILE
###################################################
## Llegim les expres d'un mateix fitxer (normalized.filtered.csv)
## ignorant la primera columna, que seran els affy IDs
exprs <- read.csv2("normalized.filtered.csv", quote = "", dec = ",")[-1]
head(exprs)
###################################################
###################################################



###################################################
###################################################
##step 2: MAP MOUSE SYMBOLS TO HUMAN GENE IDs
###################################################
## Aquesta part traudeix els rownames de Symbols a Gene IDs (tot mouse)

library(biomaRt)

head(exprs)
dim(exprs)

## Eliminem els symbols multiples, quedant-nos sols el primer de cada fila
exprs[15:20,1:5]
my.symbols <- unname(rapply(strsplit(as.character(exprs$Symbol), "//"), function(x) x[1], how="unlist"))
exprs$Symbol <- my.symbols
length(my.symbols)
exprs[15:20,1:5]

## Aqui associem els gene symbols de ratolí, als corresonents entrez gene ids de ratolí i als symbols i entrez ids humans, 
from.mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
to.mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
homologs <- getLDS(attributes = c("entrezgene", "mgi_symbol"), filters = "mgi_symbol", values = my.symbols, mart = from.mart, 
                   attributesL = c("entrezgene", "hgnc_symbol"), martL = to.mart)
dim(homologs)
head(homologs)

## rename mouse symbols column and use it merge data frames
colnames(homologs)[1:3] <- c("entrezgene.mm","Symbol","entrezgene.hs")
head(homologs)
exprs[15:20,1:5]
my.merged <- merge(exprs, homologs, by='Symbol')
dim(my.merged)
head(my.merged)

my.merged[which(my.merged$Symbol=="Cox5b"),]
## sort data frame in order to keep the shortest human entrez gene Id in case of duplicates
my.merged <- my.merged[order(my.merged$entrezgene.hs), ]
my.merged[which(my.merged$Symbol=="Cox5b"),]
## remove lines that are duplicated for all fields but the human entrez gene, keeping the first one
limcols <- dim(my.merged)[2]-3
my.merged <- my.merged[-which(duplicated(my.merged[1:limcols])), ]
my.merged[which(my.merged$Symbol=="Cox5b"),]

## remove lines having no human entrez id
length(which(is.na(my.merged$entrezgene.hs)))
my.merged <- my.merged[-which(is.na(my.merged$entrezgene.hs)),]
length(which(is.na(my.merged$entrezgene.hs)))
dim(my.merged)
head(my.merged)

## and aggregate those cases with the same human entrez ID 
## using aggregate by entrezgene.hs and mean function
head(my.merged)
length(my.merged$entrezgene.hs)
length(unique(my.merged$entrezgene.hs))
my.aggreg <- aggregate.data.frame(my.merged, by=list(my.merged$entrezgene.hs), mean)
length(my.aggreg$entrezgene.hs)

head(my.aggreg)
rownames(my.aggreg) <- my.aggreg$entrezgene.hs
head(my.aggreg)
tail(my.aggreg)

endcols <- c(dim(my.aggreg)[2]-2, dim(my.aggreg)[2]-1, dim(my.aggreg)[2])
exprs <- my.aggreg[-c(1,2,endcols)]
head(exprs)
dim(exprs)

####################################################################################
### Això seria per passar de symbols a entrez, però sols en ratolí
#ensembl.mm = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
##listFilters(ensembl.mm)
##listAttributes(ensembl.mm)[1:100,]
#my.ids <- getBM(attributes = c('mgi_symbol','entrezgene'), filters = 'mgi_symbol',
#                    values = my.symbols, mart = ensembl.mm)
#
#head(my.ids)
#length(which(is.na(my.ids$entrezgene)))
#tail(my.ids)
#length(which(!is.na(my.ids$entrezgene)))
##my.ids[which(my.ids$mgi_symbol=="Cox5b"),]
####################################################################################




###################################################
###################################################
##step 3: DEFINE COMPARIONS
###################################################

## primer creem el targets a partir dels colnames de la matriu d'expressions (=shortnames de les mostres)
type <- sapply(colnames(exprs), function(v){strsplit(v, '[.]')[[1]][1]}, USE.NAMES = FALSE)
treat <- sapply(colnames(exprs), function(v){strsplit(v, '[.]')[[1]][2]}, USE.NAMES = FALSE)
indiv <- sapply(colnames(exprs), function(v){strsplit(v, '[.]')[[1]][3]}, USE.NAMES = FALSE)
batch <- sapply(colnames(exprs), function(v){strsplit(v, '[.]')[[1]][4]}, USE.NAMES = FALSE)
group <- paste(type, treat, sep=".")
sampleID <- paste(type, treat, indiv, batch, sep = ".") # <- colnames(exprs)
targets <- data.frame(group, type, treat, indiv, batch)
rownames(targets) <- sampleID
targets

## comprovem (i canviem si cal ) etiquetes colnames per facilitar el seguent pas 
colnames(exprs) == sampleID
# colnames(exprs) <- sampleID

## Creem comparacions 
samp <- c("Veh.Cast", "Veh.DHT", "Veh.DHT", "Lep.Cast", "Lep.DHT", "Lep.DHT", "Lep.CTL", "Lep.Cast" ,"Lep.DHT")
namessamp <- c("Veh.Cast", "Veh.Dht", "Veh.Dht", "Lep.Cast", "Lep.Dht", "Lep.Dht", "Lep.Ctl", "Lep.Cast" ,"Lep.Dht")
ref <- c("Veh.CTL", "Veh.CTL", "Veh.Cast", "Lep.CTL", "Lep.CTL", "Lep.Cast", "Veh.CTL", "Veh.Cast" ,"Veh.DHT")
namesref <- c("Veh.Ctl", "Veh.Ctl", "Veh.Cast", "Lep.Ctl", "Lep.Ctl", "Lep.Cast", "Veh.Ctl", "Veh.Cast" ,"Veh.Dht")
paired <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
comparisons <- data.frame(samp, namessamp, ref, namesref, paired)
comparisons
###################################################
###################################################





###################################################
###################################################
##step 4: GENE SET ANALYSIS WITH GAGE against KEGG
###################################################
library(gage)
library(pathview)

## cut-off value for the enrichment analyses 
## (only pathways with q-value < alpha will be selected)
alpha <- 0.10  # 0.05

#comparisons <- comparisons[c(1,2,6,9,12,15), ]

stime_kegg <- system.time({
  for(i in 1:dim(comparisons)[1])  {
  #foreach(i = 1:dim(comparisons)[1]) %dopar% {
    
    setwd(workingDir)
    cond1 <- comparisons[i,"ref"]
    cond2 <- comparisons[i,"samp"]
    sample.comp <- grep(paste0(cond1,"|",cond2),colnames(exprs),value = TRUE)
    data.comp <- exprs[,sample.comp]
    dim(data.comp)
    head(data.comp)
    
    comp <- NULL
    comp <- paste0(namessamp[i],".vs.",namesref[i])
    print(paste("KEGG pathways for:", comp))
    
    if (!dir.exists(paste0("results/", comp))) dir.create(paste0("results/", comp)) 
    setwd(paste0("results/", comp))
    #getwd()
    
    # grp.idx <-  rep(NA,dim(data.comp)[2])
    # grp.idx[ grep(cond1,rownames(targets[sample.comp,]))] <- namesref[i]
    # grp.idx[ grep(cond2,rownames(targets[sample.comp,]))] <- namessamp[i]
    # grp.idx <- as.factor(grp.idx)
    # grp.idx <- factor(grp.idx,c(namesref[i],  namessamp[i]))
    
    ## joint workflow with DEseq/edgeR/limma/Cufflinks merges around here
    
    samp.idx <- grep(cond2,rownames(targets[sample.comp,]))
    ref.idx <- grep(cond1,rownames(targets[sample.comp,]))
    data(kegg.gs)


    if (!comparisons$paired[i]) {
      exprs.kegg.p <- gage(data.comp, gsets = kegg.gs, ref = ref.idx,
                          samp = samp.idx, compare = "unpaired")
      exprs.d <-  data.comp[, samp.idx] - rowMeans(data.comp[, ref.idx])
    }
    if (comparisons$paired[i]) {
      exprs.kegg.p <- gage(data.comp, gsets = kegg.gs, ref = ref.idx,
                          samp = samp.idx, compare = "paired")
      exprs.d <-  data.comp[, samp.idx] - data.comp[, ref.idx]
    }
    
    ## differential expression: log2 ratio or fold change, unpaired samples
    
    print(paste0(c("exprs.d rows=","exprs.d cols="), dim(exprs.d)))
    #head(exprs.d)
    #tail(exprs.d)
    ## write conts.d as a probe for debugging
    write.csv(file = "exprs.diffs.csv", exprs.d)
    
    #################################################
    ## select up-regulated pathways (q-value < alpha)
    kegg.up <- as.data.frame(exprs.kegg.p$greater)
    #class(kegg.up)
    #head(kegg.up)
    kegg.up.selected <- kegg.up[which(kegg.up$q.val < alpha), ]
    print(paste0("kegg.up.selected: ", dim(kegg.up.selected)[1]))
    #head(kegg.up.selected)
    #class(kegg.up.selected)
    write.csv(file = "selected.KEGG.up.csv", kegg.up.selected)
    
    ## select up-regulated pathways (q-value < alpha)
    #sel.up <- exprs.kegg.p$greater[, "q.val"] < alpha & !is.na(exprs.kegg.p$greater[,"q.val"])
    sel.up <- kegg.up$q.val < alpha
    
    ## obtain Ids of top up-regulated pathways
    path.ids.up <- rownames(exprs.kegg.p$greater)[which(sel.up)]
    path.ids.up
    path.ids2.up <- substr(path.ids.up, 1, 8)
    
    if (length(path.ids2.up)>20){limIds <- 20} else {limIds <- length(path.ids2.up)}
    ## up-regulated pathways (top 20) generated by pathview
    if(dim(exprs.d)[2]>16){
      pv.out.list.up <- sapply(path.ids2.up[1:limIds], 
                               function(pid) try(pathview(
                                                #gene.data = exprs.d, pathway.id = pid,
                                                gene.data = exprs.d, pathway.id = pid, kegg.native = FALSE, pdf.size =c(10,10),
                                                low = list(gene="green", cpd="blue"), 
                                                mid = list(gene="gray", cpd="gray"), 
                                                high = list(gene="magenta", cpd="yellow"),
                                                limit = list(gene=2, cpd=2),
                                                bins = list(gene=20, cpd=20),
                                                species = "hsa"), TRUE))
    }else{
      pv.out.list.up <- sapply(path.ids2.up[1:limIds], 
                               function(pid) try(pathview(
                                                gene.data = exprs.d, pathway.id = pid,
                                                #gene.data = exprs.d, pathway.id = pid, kegg.native=FALSE, pdf.size =c(10,10),
                                                low = list(gene="green", cpd="blue"), 
                                                mid = list(gene="gray", cpd="gray"), 
                                                high = list(gene="magenta", cpd="yellow"),
                                                limit = list(gene=2, cpd=2),
                                                bins = list(gene=20, cpd=20),
                                                species = "hsa"), TRUE))
    }
    
    #################################################
    ## select dwn-regulated pathways (q-value < alpha)
    kegg.dwn <- as.data.frame(exprs.kegg.p$less)
    #class(kegg.dwn)
    #head(kegg.dwn)
    kegg.dwn.selected <- kegg.dwn[which(kegg.dwn$q.val < alpha), ]
    print(paste0("kegg.dwn.selected: ", dim(kegg.dwn.selected)[1]))
    #head(kegg.dwn.selected)
    #class(kegg.dwn.selected)
    write.csv(file = "selected.KEGG.dwn.csv", kegg.dwn.selected)
    
    ## select dwn-regulated pathways (q-value < alpha)
    #sel.dwn <- exprs.kegg.p$less[, "q.val"] < alpha & !is.na(exprs.kegg.p$less[,"q.val"])
    sel.dwn <- kegg.dwn$q.val < alpha
    
    ## obtain Ids of top dwn-regulated pathways
    path.ids.dwn <- rownames(exprs.kegg.p$less)[which(sel.dwn)]
    path.ids.dwn
    path.ids2.dwn <- substr(path.ids.dwn, 1, 8)
    
    if (length(path.ids2.dwn)>20){limIds <- 20} else {limIds <- length(path.ids2.dwn)}
    ## dwn-regulated pathways (top 20) generated by pathview
    if(dim(exprs.d)[2]>16){
      pv.out.list.dwn <- sapply(path.ids2.dwn[1:limIds], 
                               function(pid) try(pathview(
                                 #gene.data = exprs.d, pathway.id = pid,
                                 gene.data = exprs.d, pathway.id = pid, kegg.native=FALSE, pdf.size = c(10,10),
                                 low = list(gene="green", cpd="blue"), 
                                 mid = list(gene="gray", cpd="gray"), 
                                 high = list(gene="magenta", cpd="yellow"),
                                 limit = list(gene=2, cpd=2),
                                 bins = list(gene=20, cpd=20),
                                 species = "hsa"), TRUE))
    }else{
      pv.out.list.dwn <- sapply(path.ids2.dwn[1:limIds], 
                               function(pid) try(pathview(
                                 gene.data = exprs.d, pathway.id = pid,
                                 #gene.data = exprs.d, pathway.id = pid, kegg.native=FALSE, pdf.size = c(10,10),
                                 low = list(gene="green", cpd="blue"), 
                                 mid = list(gene="gray", cpd="gray"), 
                                 high = list(gene="magenta", cpd="yellow"),
                                 limit = list(gene=2, cpd=2),
                                 bins = list(gene=20, cpd=20), 
                                 species = "hsa"), TRUE))
    }
    
    
    setwd(workingDir)
  }
})[3]
stime_kegg/60
###################################################
###################################################


###################################################
###################################################
##step 5: OVERLAPPING PATHWAYS BETWEEN COMPARISONS
###################################################

## ------------------------------------------------
##  DATA FROM selected.KEGG...csv FILES 
##  MUST BE LOADED (fins que tinguem la funcio xula feta)

library(VennDiagram)

setwd(resultsDir)

## ------------------------------------------------
## OVERLAPPING OF DOWN-REGULATED KEGG PATHWAYS (COMPARISONS FROM GROUP 1)
file1 <- read.csv("Veh.Cast.vs.Veh.Ctl/selected.KEGG.dwn.csv")
list1 <- as.character(file1$X)

file2 <- read.csv("Veh.Dht.vs.Veh.Ctl/selected.KEGG.dwn.csv")
list2 <- as.character(file2$X)

file3 <- read.csv("Veh.Dht.vs.Veh.Cast/selected.KEGG.dwn.csv")
list3 <- as.character(file3$X)

mainTitle <- "Group1.Kegg.dwn overlapping"

## Creating Venn Diagram
venn.plot <- venn.diagram(list(A = list1, 
                               B = list2, 
                               C = list3),
                          category.names = c("Veh.Cast.vs.Veh.Ctl", "Veh.Dht.vs.Veh.Ctl", "Veh.Dht.vs.Veh.Cast"),
                          fill = c("grey", "orange", "blue"),
                          alpha = 0.50,
                          resolution = 600,
                          cat.cex = 0.9,
                          main = mainTitle,
                          filename = NULL)
pdf(paste(mainTitle,"pdf",sep="."))
grid.draw(venn.plot)
dev.off()
rm(file1, file2, file3)
rm(list1, list2, list3)


## ------------------------------------------------
## OVERLAPPING OF UP-REGULATED KEGG PATHWAYS (COMPARISONS FROM GROUP 1)
file1 <- read.csv("Veh.Cast.vs.Veh.Ctl/selected.KEGG.up.csv")
list1 <- as.character(file1$X)

file2 <- read.csv("Veh.Dht.vs.Veh.Ctl/selected.KEGG.up.csv")
list2 <- as.character(file2$X)

file3 <- read.csv("Veh.Dht.vs.Veh.Cast/selected.KEGG.up.csv")
list3 <- as.character(file3$X)

mainTitle <- "Group1.Kegg.up overlapping"

## Creating Venn Diagram
venn.plot <- venn.diagram(list(A = list1, 
                               B = list2, 
                               C = list3),
                          category.names = c("Veh.Cast.vs.Veh.Ctl", "Veh.Dht.vs.Veh.Ctl", "Veh.Dht.vs.Veh.Cast"),
                          fill = c("grey", "orange", "blue"),
                          alpha = 0.50,
                          resolution = 600,
                          cat.cex = 0.9,
                          main = mainTitle,
                          filename = NULL)
pdf(paste(mainTitle,"pdf",sep="."))
grid.draw(venn.plot)
dev.off()
rm(file1, file2, file3)
rm(list1, list2, list3)


## ------------------------------------------------
## OVERLAPPING OF DOWN-REGULATED KEGG PATHWAYS (COMPARISONS FROM GROUP 2)
file1 <- read.csv("Lep.Cast.vs.Lep.Ctl/selected.KEGG.dwn.csv")
list1 <- as.character(file1$X)

file2 <- read.csv("Lep.Dht.vs.Lep.Ctl/selected.KEGG.dwn.csv")
list2 <- as.character(file2$X)

file3 <- read.csv("Lep.Dht.vs.Lep.Cast/selected.KEGG.dwn.csv")
list3 <- as.character(file3$X)

mainTitle <- "Group2.Kegg.dwn overlapping"

## Creating Venn Diagram
venn.plot <- venn.diagram(list(A = list1, 
                               B = list2, 
                               C = list3),
                          category.names = c("Lep.Cast.vs.Lep.Ctl", "Lep.Dht.vs.Lep.Ctl", "Lep.Dht.vs.Lep.Cast"),
                          fill = c("grey", "orange", "blue"),
                          alpha = 0.50,
                          resolution = 600,
                          cat.cex = 0.9,
                          main = mainTitle,
                          filename = NULL)
pdf(paste(mainTitle,"pdf",sep="."))
grid.draw(venn.plot)
dev.off()
rm(file1, file2, file3)
rm(list1, list2, list3)


## ------------------------------------------------
## OVERLAPPING OF UP-REGULATED KEGG PATHWAYS (COMPARISONS FROM GROUP 2)
file1 <- read.csv("Lep.Cast.vs.Lep.Ctl/selected.KEGG.up.csv")
list1 <- as.character(file1$X)

file2 <- read.csv("Lep.Dht.vs.Lep.Ctl/selected.KEGG.up.csv")
list2 <- as.character(file2$X)

file3 <- read.csv("Lep.Dht.vs.Lep.Cast/selected.KEGG.up.csv")
list3 <- as.character(file3$X)

mainTitle <- "Group2.Kegg.up overlapping"

## Creating Venn Diagram
venn.plot <- venn.diagram(list(A = list1, 
                               B = list2, 
                               C = list3),
                          category.names = c("Lep.Cast.vs.Lep.Ctl", "Lep.Dht.vs.Lep.Ctl", "Lep.Dht.vs.Lep.Cast"),
                          fill = c("grey", "orange", "blue"),
                          alpha = 0.50,
                          resolution = 600,
                          cat.cex = 0.9,
                          main = mainTitle,
                          filename = NULL)
pdf(paste(mainTitle,"pdf",sep="."))
grid.draw(venn.plot)
dev.off()
rm(file1, file2, file3)
rm(list1, list2, list3)


## ------------------------------------------------
## OVERLAPPING OF DOWN-REGULATED KEGG PATHWAYS (COMPARISONS FROM GROUP 3)
file1 <- read.csv("Lep.Ctl.vs.Veh.Ctl/selected.KEGG.dwn.csv")
list1 <- as.character(file1$X)

file2 <- read.csv("Lep.Cast.vs.Veh.Cast/selected.KEGG.dwn.csv")
list2 <- as.character(file2$X)

file3 <- read.csv("Lep.Dht.vs.Veh.Dht/selected.KEGG.dwn.csv")
list3 <- as.character(file3$X)

mainTitle <- "Group3.Kegg.dwn overlapping"

## Creating Venn Diagram
venn.plot <- venn.diagram(list(A = list1, 
                               B = list2, 
                               C = list3),
                          category.names = c("Lep.Ctl.vs.Veh.Ctl", "Lep.Cast.vs.Veh.Cast", "Lep.Dht.vs.Veh.Dht"),
                          fill = c("grey", "orange", "blue"),
                          alpha = 0.50,
                          resolution = 600,
                          cat.cex = 0.9,
                          main = mainTitle,
                          filename = NULL)
pdf(paste(mainTitle,"pdf",sep="."))
grid.draw(venn.plot)
dev.off()
rm(file1, file2, file3)
rm(list1, list2, list3)


## ------------------------------------------------
## OVERLAPPING OF UP-REGULATED KEGG PATHWAYS (COMPARISONS FROM GROUP 3)
file1 <- read.csv("Lep.Ctl.vs.Veh.Ctl/selected.KEGG.up.csv")
list1 <- as.character(file1$X)

file2 <- read.csv("Lep.Cast.vs.Veh.Cast/selected.KEGG.up.csv")
list2 <- as.character(file2$X)

file3 <- read.csv("Lep.Dht.vs.Veh.Dht/selected.KEGG.up.csv")
list3 <- as.character(file3$X)

mainTitle <- "Group3.Kegg.up overlapping"

## Creating Venn Diagram
venn.plot <- venn.diagram(list(A = list1, 
                               B = list2, 
                               C = list3),
                          category.names = c("Lep.Ctl.vs.Veh.Ctl", "Lep.Cast.vs.Veh.Cast", "Lep.Dht.vs.Veh.Dht"),
                          fill = c("grey", "orange", "blue"),
                          alpha = 0.50,
                          resolution = 600,
                          cat.cex = 0.9,
                          main = mainTitle,
                          filename = NULL)
pdf(paste(mainTitle,"pdf",sep="."))
grid.draw(venn.plot)
dev.off()
rm(file1, file2, file3)
rm(list1, list2, list3)


###################################################
###################################################
#### END OF THE SCRIPT
###################################################
###################################################
