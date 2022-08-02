gc()
# dev.off(dev.list()["RStudioGD"])

library(ggplot2)
library(ggpubr)
library(dplyr)
# library(mgcv)
library(tidymv)
library(patchwork) # graph organization within a figure
# library(gtsummary)
# library(grid)
# library(stringr)
# library(gridExtra)
# library(tables)
# library(grid)
# library(gridExtra)
# library(data.table)
# library(formattable)
# library(tidyr)
# library(ggseg)

#-------------------------------------------------------------------------------
# FUNCTION DEFINITION
#-------------------------------------------------------------------------------
generateScatterPlot <- function(df, phenotype, phenotypeName){
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  plot01 <- df %>%
    ggplot(aes(x=age_in_years, y=phenotype), log='x') +
    geom_point(data=df, alpha=0.5, aes(x=age_in_years, y=phenotype, color=top_scan_reason_factors)) +
    # geom_vline(xintercept=peakAge, linetype='dashed', color='#c23400') +
    scale_color_manual(values = cbbPalette, name = "Sex") +
    # scale_x_continuous(trans='log10') +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(title = phenotypeName)
  
  print(plot01)
  return(plot01)
}


#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------

# Load the data
# fn = '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-07-29_highres_nocontrast_singlescanpersubject.csv'
fn = '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_combatted_covariates_removed_plus_metadata.csv'
# t <- "Scatterplots of Precombat Phenotypes"
t <- "Scatterplots of Combatted Phenotypes (Covariate Effects Removed)"
analysisDf <- read.csv(fn)


# Print the column names
print(colnames(analysisDf))

# Identify a subset of columns
cols = c("VentricleVolume", "CorticalSurfaceArea", "TCV", "MeanCorticalThickness",
         "CortexVol", 'CerebralWhiteMatterVol', "TotalGrayVol", "eTIV", "SurfaceHoles")

# For each column in the subset, generate an age x phenotype scatterplot
scatterPlots <- c()
for (col in cols){
  print(col)
  scatterPlots[[col]] <- generateScatterPlot(analysisDf, analysisDf[ , col], col)
}

# Show the plots
patch <- wrap_plots(scatterPlots, guides='collect')
print(patch + plot_annotation(title=t))

      