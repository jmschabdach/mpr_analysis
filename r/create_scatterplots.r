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
# fn <- '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-07-29_highres_nocontrast_singlescanpersubject.csv'
# fn <- '/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
# t <- "Scatterplots of Original Phenotypes"
# fnOut <- '/Users/youngjm/Data/clip/figures/2022-09-16_trajectories_raw_phenotype_centiles.png'

# fn <- '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_combatted_plus_metadata.csv'
# t <- "Scatterplots of Combatted Phenotypes"
# fnOut <- '/Users/youngjm/Data/clip/figures/2022-09-16_trajectories_combatted_phenotypes.png'
# cols = c("TotalGrayVol", "CerebralWhiteMatterVol", "SubCortGrayVol",
#          "VentricleVolume", "CorticalSurfaceArea", "MeanCorticalThickness",
#          "TCV")

fn <- "/Users/youngjm/Data/clip/fs6_stats/combatted_phenotypes_normalized.csv"
# fnOut <- "/Users/youngjm/Data/clip/figures/2022-09-08_normalized_combatted_phenotypes.png"
# t <- "Scatterplots of Normalized Combatted Phenotypes"

# fn <- '/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_normalized.csv'
# t <- "Scatterplots of Brainchart-normalized Phenotypes"
# fnOut <- '/Users/youngjm/Data/clip/figures/2022-09-16_trajectories_original_normalized_phenotypes.png'

# cols <- c("GMVTransformed.normalised", "WMVTransformed.normalised",
#           "sGMVTransformed.normalised", "VentriclesTransformed.normalised",
#           "totalSA2Transformed.normalised", "meanCT2Transformed.normalised",
#           "TCVTransformed.normalised")

fnOut <- "/Users/youngjm/Data/clip/figures/2022-09-16_normalized_combatted_centiles.png"
t <- "Scatterplots of Centiles"

# fnOut <- "/Users/youngjm/Data/clip/figures/2022-09-16_trajectories_original_normalized_centiles.png"
# t <- "Scatterplots of Brainchart-Normalized Centiles"
# 
cols = c("GMVTransformed.q.wre", "WMVTransformed.q.wre",
         "sGMVTransformed.q.wre", "VentriclesTransformed.q.wre",
         "meanCT2Transformed.q.wre", "totalSA2Transformed.q.wre",
         "TCVTransformed.q.wre")

# fn <- "/Users/youngjm/Data/clip/images/derivatives/synthseg_2.0_phenotypes.csv"
# fnOut <- "/Users/youngjm/Data/clip/figures/2022-09-08_synthseg_global_phenotypes.png"
# t <- "Scatterplots of SynthSeg 2.0 Phenotypes"
# cols = c("GMV", "WMV",
#          "sGMV", "Ventricles",
#          # "totalSA2Transformed.q.wre", "meanCT2Transformed.q.wre",
#          "TCV")

analysisDf <- read.csv(fn)

# Synthseg output only
if ('age_in_days' %in% colnames(analysisDf)){
  names(analysisDf)[names(analysisDf) == 'age_in_days'] <- 'age_at_scan_days'
}
analysisDf <- addPrimaryScanReasonCol(analysisDf)

# Print the column names
print(colnames(analysisDf))

# Identify a subset of columns
# cols = c("VentricleVolume", "CorticalSurfaceArea", "TCV", "MeanCorticalThickness",
#          "CortexVol", 'CerebralWhiteMatterVol', "TotalGrayVol", "eTIV", "SurfaceHoles")

# For each column in the subset, generate an age x phenotype scatterplot
scatterPlots <- c()
for (col in cols){
  print(col)
  scatterPlots[[col]] <- generateScatterPlot(analysisDf, analysisDf[ , col], col)
}

# Show the plots
patch <- wrap_plots(scatterPlots, guides='collect')
png(file=fnOut,
    width=1200, height=800)
print(patch + plot_annotation(title=t))
dev.off()

# Save with figure size 
      