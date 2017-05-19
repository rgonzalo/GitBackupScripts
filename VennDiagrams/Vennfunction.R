###Aquesta funció serveix per fer diagrames de venns de dos llistes
###i també genera un .csv indicant els gens que hi ha a cada venn circle
###Queda millorar-la per tal de poder fer-ho amb més llistes
##Paràmetres
#file1 = data.frame de la primera topTable (si és amb format "BasicPipe" millor)
#name1 = nom de la comparació 1
#file2 = data.frame de la segona topTable (si és amb format "BasicPipe" millor)
#name2 = nom de la comparació 2
#FC = fold change dessitjat (el que estigui codificat a la columna "logFC")
#pval = pvalor dessitjat (el que estigui codificat a la columna "P.Value")

  vennfrom2csv<-function(file1,name1,file2,name2,FC,pval) {
    ##seleccionem els cutoffs
    ##llista1
    pass1 <- file1[which(file1$P.Value<pval & abs(file1$logFC)>FC),]
    #dim(pass1)
    ##llista2
    pass2 <- file2[which(file2$P.Value<pval & abs(file2$logFC)>FC),]
    #dim(pass2)
    
    ####mirem que no entrades repetides per Affy ID i s'ajunten
    #list1 <- as.character(sort(unique(pass1$X)))
    #length(list1)
    #list2 <- as.character(sort(unique(pass2$X)))
    #length(list2)
    
    list1 <- as.character(rownames(pass1))
    #length(list1)
    list2 <- as.character(rownames(pass2))
    #length(list2)
    
    list <- c(list1, list2)
    #length(list)
    list <- sort(unique(list))
    #length(list)
    
    ####es crea un data.frame que omplim de 0 i després de 1 si hi coexisteixen en les dues llistes
    df <- data.frame(genes = list, l1 = rep(0,length(list)), l2 = rep(0,length(list)))
    #head(df)
    
    df$l1 <- as.numeric((df$genes %in% list1)==T)
    df$l2 <- as.numeric((df$genes %in% list2)==T)
    
    ##Fem el diagrama de venn
    overlap<-calculate.overlap(x=list("list1"=list1,"list2"=list2))
    
    draw.pairwise.venn(length(overlap$a1),length(overlap$a2),length(overlap$a3),
                       category=c(name1,name2),scaled = TRUE,euler.d = TRUE, 
                       fill = c("blue", "red"),lty = "blank",cat.pos = c(190, 190))
    
    ##es grava l'arxiu de text on s'indica quin gen hi és a cada venn
    datos<-file2[which(rownames(file2) %in% df$genes),]
    datos<-datos[,-c(3:length(colnames(datos)))]
    rownames(df)<-df$genes
    datos2<-merge(datos,df,by=0)
    datos2<-datos2[,-c(3,4)]
    colnames(datos2)<-c("AffyID","Symbols",name1,name2)
    write.csv(datos2, file=file.path(resultsDir,paste("mult.comp.",name1,"_",name2,".csv",sep = "")),sep=";")
  }