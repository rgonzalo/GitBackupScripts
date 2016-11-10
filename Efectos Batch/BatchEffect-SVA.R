################################
########      SVA       ########
################################
### Identificar variables
### latents 
################################
## http://www.bioconductor.org/packages/release/bioc/manuals/sva/man/sva.pdf
################################

if (!require(sva)) install.packages("sva") 


load("dadesexemple/normalizedData2013.Rda") ## llegim ExpresionSet de dades normalitzades
targets2013 <- read.csv2 ("dadesexemple/targets2013.csv",sep="\t",row.names=1) ## llegim targets
pData(my.norm) <- targets2013  ## Posem el nostre targets

false <- which(rownames(pData(my.norm)) %in% colnames(exprs(my.norm))==FALSE) # optatiu, només en el cas d'eliminar alguna mostra
pheno <- pData(my.norm)[-c(false),] ## si no s'elimina cap mostra treure []
edata <- exprs(my.norm)[,rownames(pheno)]

mod <- model.matrix ( ~ as.factor(Grupo) , data = pheno) ## Totes les variables
mod0 <- model.matrix( ~ 1, data = pheno) ## Variables que poden afectar sense la nostre variable d'interes


##### Estimem diferents Batch
## First it identifies the number of latent factors that need to be estimate
n.sv = num.sv(edata,mod,method="leek")
n.sv                                  ## numero de variables latents
## estimate surrogate variables
svobj = sva(edata,mod,mod0,n.sv=n.sv)


##### Adjusting for surrogate variables using the f.pvalue function
pValues = f.pvalue(edata,mod,mod0)
(sum(pValues < 0.01)/length(pValues))*100 ## percentatge p.valors significatius
qValues = p.adjust(pValues,method="BH")
(sum(qValues < 0.01)/length(qValues))*100 ## percentatge p.val.adj significatius

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)    ## p.value surrogate variables
qValuesSv = p.adjust(pValuesSv,method="BH") ## p.val.adj surrogate variables


##### Adjusting for surrogate variables using the limma package
fit = lmFit(edata,modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)


##### Applying the ComBat function to adjust for known batches
batch = c(pheno$Gender,pheno$Batch)
modcombat = model.matrix(~1, data=pheno)

## matriu d'expresions ajustada per variables batch
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, numCovs=NULL, par.prior=TRUE,prior.plots=TRUE)

pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

