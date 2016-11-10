############################################################################
##
## PARAMETERS SOURCE FILE TO BE USED BY THE SCRIPT generateHeatMap.R
## Should be included in generateHeatMap.R with source("heatMapParams.R")
##
## Created: 2013-04-29 by Ferran Brianso (ferran.brianso@vhir.org)
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
## PARAMETERS REGARDING INPUT FILE WITH EXPRESSION DATA IN CSV FORMAT
############################################################################
csv.file <- "inputExprValues.csv" # edit input csv file name
dec.sep <- "," # set current decimal separator in the csv file
col.sep <- ";" # set current column (field) separator in the csv file
ids.col <- 1 # set column index where to find the feature labels (e.g. GeneSymbols)
############################################################################



############################################################################
## PARAMETERS REGARDING INPUT LIST OF IDS TO BE INCLUDED IN THE HEAT MAP
############################################################################
sublist.file <- "subsetList.txt" # edit list of row IDs to be included in the heat map (should match some IDs in input csv file)
                                 # or leave as sublist.file <- "" if no subsetting desired
label.lim <- 6 # string position at which to trim feature labels being included in the plot (0 if full label allowed)
############################################################################



############################################################################
## PARAMETERS REGARDING SAMPLES WITH CORRESPONDING NAMES, LABELS AND COLORS
############################################################################
samples2clust <- c(2:6,17:21,23:27,28:31) # edit indexes of columns corresponding with samples to be included in the heat map
                                          # (use , to separate and : to include all in range)          
replace.labels <- TRUE # set to TRUE to use col.names instead of the csv file header for samples
col.names <- c("CTL.A1","CTL.A2C","CTL.A5","CTL.A7B","CTL.A3B", # edit labels to use in the heat map (instead of column headers)
               "DM.B3B","DM.B2","DM.E8B","DM.B9","DM.B6",
               "PM.E1","PM.E6","PM.B5","PM.E7C","PM.E5",
               "IBM.F7B","IBM.F6B","IBM.F4B","IBM.F2")
colors4samples <- c(rep("green", 5),  # edit color and number of times to repeat it
                    rep("red", 5),    # idem
                    rep("magenta", 5),# idem
                    rep("violet", 4)) # idem
############################################################################



############################################################################
## PARAMETERS REGARDING THE HEAT MAP GENERATION
############################################################################
dendro.rows <- TRUE # make Y axis dendrogram to order features by expression profile (TRUE by default)
dendro.cols <- FALSE # make X axis dendrogram to order samples by expression profile (FALSE by default)
title.text <- "Heat map for DM, PM and IBM vs. CTRL \n samples (genes with p-value < 0.01)" # edit heat map title
row.cex <- 0.6 #cex numeric value for row feature labels (0.75 by default)
col.cex <- 0.75 #cex numeric value for column sample labels (0.75 by default)
color.palette <- colorpanel(32, "green", "black", "red") # choose 3 colors (for negative, 0 and positive values)
num.clusters <- 0 # a priori number of groups for feature clustering (if 0, param colors4samples is used to compute it)
############################################################################



############################################################################
## PARAMETERS REGARDING DETAILS OF THE HEAT MAP OUTPUT FILES
############################################################################
out.name <- "heatmap" # edit heat map file name without extension
out.type <- "tif" # set heat map desired extension ("eps", "tif", "bmp" or "png" allowed) in addition to default pdf (put "pdf" if only pdf format desired)
out.dpi <- 300 # set the dpi (dots per inch) density wanted for the output heat map. Set to 0 if no dpi requirements
############################################################################




