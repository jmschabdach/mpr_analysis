gc()
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")
# dev.off(dev.list()["RStudioGD"])

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidymv)
library(patchwork) # graph organization within a figure

library(mgcv)
library(gtsummary)
library(grid)
library(stringr)
library(gridExtra)
library(reshape2)
library(tables)
library(grid)
library(gridExtra)
library(data.table)
library(formattable)
library(tidyr)
library(ggseg)


#-------------------------------------------------------------------------------
# FUNCTION DEFINITION
#-------------------------------------------------------------------------------
generateScatterPlot <- function(df, phenotype, phenotypeName){
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  plot01 <- df %>%
    ggplot(aes(x=age_at_scan_days, y=phenotype), log='x') +
    geom_point(data=df, alpha=0.5, aes(x=age_at_scan_days, y=phenotype, color=top_scan_reason_factors)) +
    scale_color_manual(values = cbbPalette, name = "Reason for Scan") +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(title = phenotypeName)
  
  print(plot01)
  return(plot01)
}


#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------

# Load the data
# phenotypes = c("GMV", "WMV", "sGMV", "Ventricles", "TCV")
# pheno = "GMV"
# cols = c("TotalGrayVol", "TotalGrayVol", paste0(pheno, 'Transformed.normalised'), paste0(pheno, 'Transformed.normalised'))
# pheno = "WMV"
# cols = c("CerebralWhiteMatterVol", "CerebralWhiteMatterVol", paste0(pheno, 'Transformed.normalised'), paste0(pheno, 'Transformed.normalised'))
# pheno = "sGMV"
# cols = c("SubCortGrayVol", "SubCortGrayVol", paste0(pheno, 'Transformed.normalised'), paste0(pheno, 'Transformed.normalised'))
# pheno = "Ventricles"
# cols = c("VentricleVolume", "VentricleVolume", paste0(pheno, 'Transformed.normalised'), paste0(pheno, 'Transformed.normalised'))
# pheno = "CT"
# cols = c("MeanCorticalThickness", "MeanCorticalThickness", 'meanCT2Transformed.normalised', 'meanCT2Transformed.normalised')
pheno = "SA"
cols = c("CorticalSurfaceArea", "CorticalSurfaceArea", 'totalSA2Transformed.normalised', 'totalSA2Transformed.normalised')
# pheno = "TCV"
# cols = c("TCV", "TCV", paste0(pheno, 'Transformed.normalised'), paste0(pheno, 'Transformed.normalised'))


fnOrig <- '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-07-29_highres_nocontrast_singlescanpersubject.csv'
tOrig <- paste0("Scatterplot of \nPrecombat Phenotypes for ", pheno)
fnCombat <- '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_combatted_plus_metadata.csv'
tCombat <- paste0("Scatterplot of \nCombatted Phenotypes for ", pheno)
fnCNormalized <- "/Users/youngjm/Data/clip/fs6_stats/combatted_phenotypes_normalized.csv"
tCNormalized <- paste0("Scatterplot of \nCombat+Brainchart Normalized Phenotypes for ", pheno)
fnNormalized <- "/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_normalized.csv"
tNormalized <- paste0("Scatterplot of \nBrainchart Normalized Phenotypes for ", pheno)
# fnSynthSeg <- "/Users/youngjm/Data/clip/images/derivatives/synthseg_2.0_phenotypes.csv"
# tSynthSeg <- paste0("Scatterplot of SynthSeg 2.0 Phenotypes for ", pheno)

fnOut <- paste0("/Users/youngjm/Data/clip/figures/2022-09-19_phenotype_scatterplots_",pheno,".png")
scatterPlots <- c()

# For each column in the subset, generate an age x phenotype scatterplot
# Original phenotype
analysisDf <- read.csv(fnOrig)
analysisDf <- addPrimaryScanReasonCol(analysisDf)
scatterPlots[[1]] <- generateScatterPlot(analysisDf, analysisDf[ , cols[1]], tOrig)

# Combatted phenotype - this one seems off, is this the correct file name?
analysisDf <- read.csv(fnCombat)
analysisDf <- addPrimaryScanReasonCol(analysisDf)
scatterPlots[[2]] <- generateScatterPlot(analysisDf, analysisDf[ , cols[2]], tCombat)

# Normalized phenotype
analysisDf <- read.csv(fnCNormalized)
analysisDf <- addPrimaryScanReasonCol(analysisDf)
scatterPlots[[3]] <- generateScatterPlot(analysisDf, analysisDf[ , cols[3]], tCNormalized)

# Normalized phenotype
analysisDf <- read.csv(fnNormalized)
analysisDf <- addPrimaryScanReasonCol(analysisDf)
scatterPlots[[4]] <- generateScatterPlot(analysisDf, analysisDf[ , cols[4]], tNormalized)

# # SynthSeg phenotype
# analysisDf <- read.csv(fnSynthSeg)
# # Drop any rows where neurofibromatosis is in the scan_reason_categories column
# analysisDf <- masterDf[!grepl("neurofibromatosis", masterDf$scan_reason_primary), ]
# analysisDf <- analysisDf[!grepl("neurofibromatosis", analysisDf$scan_reason_categories), ]
# # Need to convert missing values in "confirm_neurofibromatosis" to FALSE
# analysisDf <- analysisDf %>%
#   mutate(confirm_neurofibromatosis = case_when(
#     grepl('1', confirm_neurofibromatosis, fixed=TRUE) ~ TRUE,
#     grepl('0', confirm_neurofibromatosis, fixed=TRUE) ~ FALSE
#   ))
# analysisDf <- analysisDf[(analysisDf$confirm_neurofibromatosis == FALSE), ]
# analysisDf <- addPrimaryScanReasonCol(analysisDf)
# # analysisDf$top_scan_reason_factors <- as.factor(analysisDf$top_scan_reason_factors)
# scatterPlots[[4]] <- generateScatterPlot(analysisDf, analysisDf[ , cols[4]], tSynthSeg)



# Show the plots
patch <- wrap_plots(scatterPlots, ncol=length(scatterPlots), guides='collect')
png(file=fnOut,
    width=1500, height=400)
print(patch + plot_annotation(title=pheno))
dev.off()

# Save with figure size 
