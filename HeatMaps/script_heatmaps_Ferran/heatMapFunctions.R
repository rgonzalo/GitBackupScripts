############################################################################
##
## PARAMETERS SOURCE FILE TO BE USED BY THE SCRIPT generateHeatMap.R
## Should be included in generateHeatMap.R with source("heatMapFunctions.R")
##
## Created: 2011-04-19 by Alex Sanchez (alex.sanchez@vhir.org)
## Modified: 2013-04-30 by Ferran Brianso (ferran.brianso@vhir.org)
## More info: http://ueb.vhir.org (ueb@vhir.org)
##
## Dual Licensed under the following copyright licenses:
## a) Creative Commons Attribution-ShareAlike 3.0 License : cc-by-sa 3.0
## b) GNU/LGPL : http://www.gnu.org/copyleft/lesser.html
##
## Acknowledgements: 
## * to all UEB staff for comments and feedback
############################################################################
############################################################################



# ------------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------------
library(gplots)

# ------------------------------------------------------------------------
# function distanceMat
# ------------------------------------------------------------------------
distanceMat <- function (x, distfun){
if (is.null(distfun)){
  stop('A distance function must be provided')
  }else{
    if (distfun == 'cor'){
      dMat <- as.dist(1 - cor(t(x)))
    }else{
      dMat <- dist(t(x), method=distfun)
    }
  }
  return(dMat)
}

# ------------------------------------------------------------------------
# function clusterAnalysis
# ------------------------------------------------------------------------
clusterAnalysis <- function(expres,
                            genes,
                            samples,
                            sampleNames,
                            comparisonName,
                            anotPackage,
                            outputDir,
                            toPDF=TRUE,
                            numClusters,
                            rowDistance='cor',
                            colDistance='euclidean',
                            RowVals = TRUE,
                            ColVals = TRUE,
                            colorsSet,
                            colsForGroups,
                            escala,
                            densityInfo = "none",
                            cexForColumns,
                            cexForRows,
                            plotProfiles = FALSE,
                            Title = "",
                            csvType = NULL)
{
  exprs2Cluster <- expres[genes, samples]
  if (!( is.null(sampleNames)))
    colnames(exprs2Cluster) <- sampleNames

  ### 1st step: heat map representation

  # dendrogram settings
  dendro <- 'both'

  if (RowVals && (!is.numeric(RowVals)))
  {
      distMat <- distanceMat(exprs2Cluster, rowDistance)
      clustRow <- hclust(distMat, method = "average")
      dendroRow <- as.dendrogram(clustRow)
  }else{
      dendroRow <- RowVals
  }

  if (ColVals && (!is.numeric(ColVals)))
  {
    if(!is.null(colDistance)){
      distMat <- distanceMat(exprs2Cluster, colDistance)
      clustCol <- hclust(distMat, method = "average")
      dendroCol <-  as.dendrogram(clustCol)
      }else{
          dendroCol <- ColVals
      }
  }else{
    if(is.numeric(ColVals))
      {
        exprs2Cluster <- exprs2Cluster[, ColVals]
        colsForGroups <- colsForGroups[ColVals]
        ColVals <- FALSE
        dendroCol <- FALSE
      }else{
        dendroCol <- ColVals
      }
  }

  if (ColVals && RowVals)
  {
    dendro <- 'both'
  }else{
    if (ColVals && (!RowVals))
    {
      dendro <- 'column'
    }else{
      if ((!ColVals) && RowVals)
      {
        dendro <- 'row'
      }else{
        dendro <- 'none'
      }
    }
  }

  mainTitle <- ifelse (Title == "", comparisonName, Title )

  if(toPDF)
  {
    heatMapFName <- paste("HeatMap", comparisonName, "pdf", sep=".")
    pdf(file.path(outputDir, heatMapFName))
  }

  hm <- heatmap.2(exprs2Cluster,
                  col = colorsSet,
                  ColSideColors = as.character(colsForGroups),
                  scale = escala,
                  Rowv = dendroRow,
                  Colv = dendroCol,
                  dendrogram = dendro,
                  key = TRUE,
                  symkey = FALSE,
                  density = densityInfo,
                  trace = "none",
                  cexCol = cexForColumns,
                  cexRow = cexForRows,
                  main = mainTitle)

  if(toPDF)
  {
    dev.off()
  }

  ### 2nd step: if some clusters given (!is.null(numClusters))
  ###           do the average profile for each cluster
  ###           write file mapping genes to clusters

  if  (RowVals)
  {
    cutN <- cutree(clustRow, numClusters)

    if (ColVals)
    {
      names.ord <- (clustCol$labels[clustCol$order])
    }else{
      names.ord <- 1:ncol(exprs2Cluster)
    }

    if(plotProfiles)
    {
    if(toPDF) {
      pdf(file.path(outputDir, plotClustFName))
      plotClustFName <- paste("ProfilePlots", comparisonName, "pdf", sep = ".")
      }
         for(i in 1:numClusters)
         {
            plotClust(exprs2Cluster[, names.ord], cutN, i, scal = T)
          }
    if(toPDF) dev.off()
    }

  ### 3rd step: extract gene IDs (gene symbols) and write all into the csv

  gNames <- rownames(exprs2Cluster[hm$rowInd, ])
  if (!is.null(anotPackage)){
    my_SYMBOL_env <- eval(parse(text = paste(anotPackage, "SYMBOL",sep = "")))
    mySymbols <- unlist(mget(gNames, my_SYMBOL_env))
    my.genes <- data.frame(symbol = mySymbols, ID = gNames, GrupID = cutN[gNames], exprs2Cluster[gNames, names.ord])
  }else{
    my.genes <- data.frame(ID = gNames, GrupID = cutN[gNames], exprs2Cluster[gNames, names.ord])
  }

  genesInClustersFName <- paste("genesInClusters", comparisonName, sep = ".")
  write.table(my.genes, file = file.path(outputDir, genesInClustersFName), sep="\t", row.names = FALSE, quote = FALSE)
  }

  return(hm)
}


# ------------------------------------------------------------------------
# function new.clustPar: Creates a clustPar object with given parameters
# ------------------------------------------------------------------------
new.clustPar <- function (expres,
                       expresFileName  = NULL,
                       genes2cluster = NULL,
                       geneListFName  = NULL,
                       samples2cluster = NULL,
                       samplesListFName  = NULL,
                       sampleNames = NULL,
                       comparisonName = "",
                       anotPackage = NULL,
                       outputDir =  ".",
                       toPDF = FALSE,
                       numClusters = 2,
                       rowDistance = 'euclidean',
                       colDistance = 'cor',
                       RowVals = TRUE,
                       ColVals = TRUE,
                       escala = "row",
                       colorsSet =  NULL,
                       densityInfo = "density",
                       colsForGroups =NULL,
                       cexForColumns = 0.8,
                       cexForRows = 0.8,
                       Title = NULL)
{
    clustPar <- list ( expres = expres ,
                       expresFileName  = expresFileName,
                       genes2cluster = genes2cluster,
                       geneListFName  = geneListFName,
                       samples2cluster = samples2cluster,
                       samplesListFName  = samplesListFName,
                       sampleNames = sampleNames,
                       comparisonName = comparisonName,
                       anotPackage = anotPackage,
                       outputDir =  outputDir,
                       toPDF = toPDF,
                       numClusters = numClusters,
                       rowDistance = rowDistance,
                       colDistance = colDistance,
                       RowVals = RowVals,
                       ColVals = ColVals,
                       escala = escala,
                       colorsSet =  colorsSet,
                       densityInfo = densityInfo,
                       colsForGroups =colsForGroups,
                       cexForColumns = cexForColumns,
                       cexForRows = cexForRows,
                       Title = Title)
    return(clustPar)
}


# ------------------------------------------------------------------------
# function doClusterAnalysis : Runs 'clusterAnalysis' with clustPar contents
# ------------------------------------------------------------------------
doClusterAnalysis <- function(clustPar)
{

  p <- clustPar

  if (!is.null(p$expresFileName)){
      expres <- loadFromFile (file.path(p$outputDir, p$expresFileName))
  }else{
      if (!is.null(p$expres)) {
	       expres <- eval(parse(text = p$expres))
      }else{
        	stop("Error, file name or data not defined")
      }
  }

  if (!is.null(p$geneListFName)){
      genes2cluster <- loadFromFile (file.path(p$outputDir, p$geneListFName))
  }else{
      if (is.null(p$genes2cluster)) {
         genes2cluster <-rownames(expres)
      }else{
        genes2cluster <- p$genes2cluster
      }
  }

  if (is.null(p$samples2cluster)) {
         samples2cluster <-1:ncol(expres)
      }else{
        samples2cluster <- p$samples2cluster
      }
      
  clust <- clusterAnalysis(expres = expres,
                           genes = genes2cluster,
                           samples = samples2cluster,
                           sampleNames = p$sampleNames,
                           comparisonName = p$comparisonName,
                           anotPackage = p$anotPackage,
                           outputDir = p$outputDir,
                           toPDF = p$toPDF, 
                           numClusters = p$numClusters,
                           rowDistance = p$rowDistance,
                           colDistance = p$colDistance,
                           RowVals = p$RowVals,
                           ColVals = p$ColVals,
                           escala = p$escala,
                           colorsSet = p$colorsSet,
                           densityInfo = p$densityInfo,
                           colsForGroups = p$colsForGroups,
                           cexForColumns = p$cexForColumns,
                           cexForRows = p$cexForRows,
                           Title = p$Title,
                           csvType = p$csvType)

  return(clust)
}


