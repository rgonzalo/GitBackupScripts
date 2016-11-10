
genesFromTopTable <- function (aTopTable,
                               filename=NULL, 
                               entrezOnly = TRUE, 
                               uniqueIds=TRUE, 
                               adjOrrawP = "adj",
                               Pcutoff = 0.05,
                               FCcutoff = 1,
                               id2Select = "Entrezs", # c("Entrezs","SYMBOL", "other")
                               cols2Select = 2){
  
  if (! is.null(filename)){
    topTab <- read.csv(file.path(listsDir, filename), head=TRUE, sep=";", dec=",", row.names=1)
  }else{
    topTab=aTopTable
  }
  dim(topTab)
  colnames(topTab)
  if (entrezOnly) {
    selectedEntrez <- !is.na(topTab[,"EntrezsA"])
    topTab <- topTab[selectedEntrez,]
  }
  dim(topTab)
  if (Pcutoff < 1){
    if (adjOrrawP=="adj"){
      selectedP <- topTab[,"adj.P.Val"] < Pcutoff
    }else{
      selectedP <- topTab[,"P.Value"] < Pcutoff
    }
    topTab<- topTab[selectedP, ]
  }
  dim(topTab)
  if (FCcutoff > 0){
    selectedFC <-(abs(topTab[,"logFC"]) > FCcutoff)
    topTab<- topTab[selectedFC, ]
    dim(topTab)
  }
  
  if (id2Select=="SYMBOL"){
    geneList <- topTab[,"SymbolsA"]
  }else{
    if (id2Select=="Entrez"){
      geneList <- topTab[,"EntrezsA"]
    }else{
      if(cols2Select> 0){
        geneList <- topTab[,cols2Select]
      }else{
        geneList <-rownames(topTab)
      }}}
  # length(geneList)
  if(uniqueIds) geneList <- unique(geneList)
  return(geneList)
}


extractInfo <- function (x, compName, pattern, outDir, adjOrraw="adj", pCutOff=0.01, fcCutoff=1){
  geneList <- genesFromTopTable (x, entrezOnly = TRUE, uniqueIds=TRUE, #uniqueEntrezs, 
                                 adjOrrawP = adjOrraw, Pcutoff = pCutOff, FCcutoff = fcCutoff, 
                                 id2Select = "EntrezsA" , cols2Select =2)
  symbolsList <- as.character(genesFromTopTable (x, entrezOnly = TRUE, uniqueIds=TRUE, #uniqueEntrezs, 
                                                 adjOrrawP = adjOrraw, Pcutoff = pCutOff, FCcutoff = fcCutoff, 
                                                 id2Select = "SymbolsA" , cols2Select =1))
  universeList <- genesFromTopTable (x, entrezOnly = TRUE, uniqueIds=TRUE, #uniqueEntrezs, 
                                     adjOrrawP = "raw", Pcutoff = 1, FCcutoff = 0, 
                                     id2Select = "EntrezsA" , cols2Select =2)
  geneList <- as.character(geneList); symbolsList <- as.character(symbolsList); universeList <- as.character(universeList);
  
  write.table(geneList, file=file.path(outDir, paste(compName,"EntrezSelected.csv", sep=".")), 
              row.names=FALSE, quote=FALSE, col.names=FALSE)
  write.table(symbolsList, file=file.path(outDir, paste(compName,"SymbolsSelected.csv", sep=".")), 
              row.names=FALSE, quote=FALSE, col.names=FALSE)
  write.table(universeList, file=file.path(outDir, paste(compName,"EntrezAll.csv", sep=".")), 
              row.names=FALSE, quote=FALSE, col.names=FALSE)
  
  columns2select <- grep(pattern, colnames(x))
  write.csv(x[,columns2select], file=file.path(outDir, paste(compName,"Expres.csv", sep=".")), 
            row.names=TRUE, quote=FALSE)
  xbyEnt<- aggregate(x[,columns2select],by= list(x$EntrezsA), mean); 
  rownames(xbyEnt)<-xbyEnt[,1];xbyEnt<-xbyEnt[,-1]
  write.csv(xbyEnt,
            file=file.path(outDir, paste(compName,"ExpresByEntrez.csv", sep=".")), 
            row.names=TRUE, quote=FALSE)
  return(list(geneList, universeList))
}

