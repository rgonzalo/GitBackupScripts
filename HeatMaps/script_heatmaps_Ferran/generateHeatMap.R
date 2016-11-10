# usage from R console: source("generateHeatMap.R")
# Created: 2013-04-29 by Ferran Brianso (ferran.brianso@vhir.org)
# Modified: 2013-04-30 by Ferran Brianso (ferran.brianso@vhir.org)
############################################################################
# SCRIPT: generateHeatMap.R
# SCOPE:  Script that creates heat maps from a single expression .csv file
#         Allows to use only a subset of genes and samples instead of the full file
# Author: Ferran Brianso (2013)
#         ferran.brianso@vhir.org
#         http://ueb.vhir.org
############################################################################
##
## Dual Licensed under the following copyright licenses:
## a) Creative Commons Attribution-ShareAlike 3.0 License : cc-by-sa 3.0
## b) GNU/LGPL : http://www.gnu.org/copyleft/lesser.html
##
##
## Acknowledgements: 
## * to all UEB staff for comments and feedback.
############################################################################


# ------------------------------------------------------------------------
# Step 0. Hide warning messages
# ------------------------------------------------------------------------
options(warn=-1)


# ------------------------------------------------------------------------
# Step 1. Load packages and specific functions
# ------------------------------------------------------------------------
source("heatMapFunctions.R") # requires R package gplots


# ------------------------------------------------------------------------
# Step 2. Load parameters from previously edited file
# ------------------------------------------------------------------------
source("heatMapParams.R") # file with all required parameters


# ------------------------------------------------------------------------
# Step 3. Read input file with expressions in a csv format
# ------------------------------------------------------------------------
x <- read.csv(file = csv.file, header = TRUE, dec = dec.sep, sep = col.sep, row.names=1)
cat(paste("\nInput", csv.file, "loaded\n", sep=" "))
cat(paste(" data frame with", dim(x)[1], "rows and", dim(x)[2], "columns\n", sep=" "))

# ------------------------------------------------------------------------
# Step 4. Subset expression data by list of given IDs and desired samples
# ------------------------------------------------------------------------
if (sublist.file != ""){
  sublist <- as.character(read.csv(sublist.file, header=F)[,1]) # list of row IDs (rownames in the input.matrix) to be included in the heat map
  cat(paste("\nInput", sublist.file, "loaded\n", sep=" "))
  cat(paste("  list with", length(sublist), "elements\n", sep=" "))
  rows2clust <- match(sublist,rownames(x))
  cat(paste("  from which", length(rows2clust), "match with data frame row names\n", sep=" "))
  cat(paste("\nSubsetting input data frame by", length(rows2clust), "rows and", length(samples2clust), "columns\n", sep=" "))
  subx <- x[rows2clust,c(ids.col,samples2clust)]
}else{
  subx <- x
}


# ------------------------------------------------------------------------
# Step 5. Put GeneSymbols+AffyIDs as row names and remove GeneSymbols column
# ------------------------------------------------------------------------
if (label.lim == 0){
  rownames(subx) <- paste(subx[,ids.col],rownames(subx),sep=".")
}else{
  rownames(subx) <- paste(substr(subx[,ids.col],1,label.lim),rownames(subx),sep=".")
}
subx <- subx[,-ids.col]


# ------------------------------------------------------------------------
# Step 6. Scale all values
# ------------------------------------------------------------------------
scaledsubx  <- scale(subx)


# ------------------------------------------------------------------------
# Step 7. Create the clustPar object with all parameters
# ------------------------------------------------------------------------
if (replace.labels == FALSE){
  col.names <- colnames(subx)
}

if (num.clusters == 0){
  num.clusters <- length(unique(colors4samples))
}

clustPar <- new.clustPar(expres = "scaledsubx" ,
                       expresFileName  = NULL,
                       genes2cluster = NULL,
                       geneListFName  = NULL,
                       samples2cluster = NULL,
                       samplesListFName  = NULL,
                       sampleNames = col.names,
                       comparisonName = "csv",
                       anotPackage = NULL,
                       outputDir =  ".",
                       toPDF = FALSE,
                       numClusters = num.clusters,
                       rowDistance = 'cor',
                       colDistance = 'euclidean',
                       RowVals = dendro.rows,
                       ColVals = dendro.cols,
                       escala = "row",
                       colorsSet = color.palette,
                       densityInfo = "density",
                       colsForGroups = colors4samples,
                       cexForColumns = col.cex, 
                       cexForRows = row.cex,
                       Title = title.text
   )


# ------------------------------------------------------------------------
# Step 8. Generate the heat map in pdf and convert to the other format
# ------------------------------------------------------------------------
name.pdf <- paste(out.name, "pdf", sep=".")
pdf(name.pdf)
hm <- doClusterAnalysis(clustPar)
dev.off()
cat(paste("\nOuput", name.pdf, "created\n", sep=" "))

if (out.type != "pdf"){
  name.new <- paste(out.name, out.type, sep=".") 
  cmd <- paste ("convert -density", out.dpi, name.pdf, name.new, sep=" ")
  system(cmd)
  cat(paste("Created", name.new, "from", name.pdf, "with ImageMagick tool\n", sep=" "))
}
  cat("\n")

# ------------------------------------------------------------------------
# END OF SCRIPT
############################################################################


