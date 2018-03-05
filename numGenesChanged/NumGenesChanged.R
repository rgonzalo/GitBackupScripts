##### topTab: Toptable en format '.csv'

####1er AMB LA FUNCIÓ "genesSelectable" ES GENEREN UNS VECTORS.
####2on FAS UN DATA.FRAME
####3er HO GUARDES EN UN .CSV

##### adj0: p.valor adjustat pel qual filtrar (de normal, 0.01)
##### adj1: segon p.valor adjustat pel qual filtrar (de normal, 0.05)
##### adj2: tercer p.valor adjustat pel qual filtrar (de normal, 0.25)
##### P1: p.valor pel qual filtrar (de normal, 0.01)
##### P2: segon p.valor pel qual filtrar (de normal, 0.05)
################################################################################

genesSelectable <- function (topTab, adj0, adj1, adj2, P1, P2,FC=0)
{
  upBelowB <- sum(topTab$B > 0  & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowB <- sum(topTab$B > 0  & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowAdj0 <- sum(topTab$adj.P.Val < adj0 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowAdj0 <- sum(topTab$adj.P.Val < adj0 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowAdj1 <- sum(topTab$adj.P.Val < adj1 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowAdj1 <- sum(topTab$adj.P.Val < adj1 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowAdj2 <- sum(topTab$adj.P.Val < adj2 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowAdj2 <- sum(topTab$adj.P.Val < adj2 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowP1 <- sum(topTab$P.Value < P1 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowP1 <- sum (topTab$P.Value < P1 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowP2 <- sum(topTab$P.Value < P2 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowP2 <- sum(topTab$P.Value < P2 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  return(c(upReg_B = upBelowB,downReg_B = downBelowB,
           upRegAdj0.01 = upBelowAdj0, downRegAdj0.01 = downBelowAdj0,
           upRegAdj0.05 = upBelowAdj1, downRegAdj0.05 = downBelowAdj1,
           upRegAdj0.25 = upBelowAdj2, downRegAdj0.25 = downBelowAdj2,
           upRegP0.01 = upBelowP1, downRegP0.01 = downBelowP1,
           upRegP0.05 = upBelowP2, downRegP0.05 = downBelowP2))   
}

####DESPRÉS S'HA DE FER EL DATA.FRAME AMB CADA VECTOR CREAT I AIXÒ ÉS EL
####GUARDES EN FORMAT .CSV
#compar1<-genesSelectable(topTab.compar1,0.01,0.05,0.25,0.01,0.05)
#numGenesChanged<-cbind.data.frame(compar1,compar2,compar3,4)
#write.csv2(numGenesChanged,file.path(resultsDir,"numGenesChanged.csv"),sep=";")