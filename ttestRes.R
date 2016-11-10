
############################
## Miriam Mota Foix
## 2016.09.07
############################

############################
## ttestRes
############################
###
### varCat = nom (de tipus caracter) de la variable categorica (binaria)
### varCons = vector amb el nom de les variables continues que es volen testar
### data = nom del data.frame on es troben les variables
### tex = TRUE / FALSE depenent si es vol el resultat mitjançant xtable
### title = en el cas de que tex sigui TRUE, caption de la taula 
###
### Generar taula amb els resultats d'un t.test, incloent mitjanes i p-valor
###
############################

ttestRes <- function(varCat, varCons, data, title = "", tex=FALSE){
res <- list()
for (i in 1:length(varCons)){
tt <- t.test( data[,varCons[i]] ~ data[ ,varCat], data = data)
res[[i]] <- c(tt$estimate[1],tt$estimate[2],tt$p.value)
}
df_sl <- data.frame(matrix(unlist(res), nrow = length(varCons), byrow = T),stringsAsFactors = FALSE)
colnames(df_sl) <- c(names(tt$estimate)[1],names(tt$estimate)[2], "P.valor")
if(!tex) return(df_sl)
if(tex) return(xtable(df_sl, caption = title))
}


############################
## Exemple
############################
fact <- as.factor(c("a","a","a","b","b","b"))
cont <- rnorm(6)
cont1 <- rnorm(6)
cont2 <- rnorm(6)
dat <- data.frame(fact,cont, cont1,cont2)

varCon <- c("cont","cont1","cont2")

ttestRes(varCat = "fact", varCons= varCon, data=dat, tex=T)
