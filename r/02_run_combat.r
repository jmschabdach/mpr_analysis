library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgcv)
library(tidymv)
library(patchwork) # graph organization within a figure
library(gtsummary)
library(grid)
# library(harrypotter)
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
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")


## Step 3: ComBat the data/Load ComBat data ------------------------------------

fnBase <- '/Users/youngjm/Data/slip/fs6_stats/'

##
# FreeSurfer ComBat

fn = '/Users/youngjm/Data/slip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
fnOut <- '/Users/youngjm/Data/slip/fs6_stats/06_combatted_fs_plus_metadata.csv'

analysisDf <- read.csv(fn)

prepForCombat(analysisDf, "/Users/youngjm/Data/slip/fs6_stats/04_toCombat_fs")

combatCommand <- paste0("python /Users/youngjm/Projects/mpr_analysis/runNeuroHarmonize.py",
                        " -p ", fnBase, "04_toCombat_fs_phenotypes.csv",
                        " -c ", fnBase, "04_toCombat_fs_covariates.csv",
                        " -o ", fnBase, "05_fs_postCombat.csv")

# Run ComBat
system(combatCommand)

# Load the combatted dataframe
combattedDf <- loadCombattedData(analysisDf, '/Users/youngjm/Data/slip/fs6_stats/05_fs_postCombat.csv')

# Keep only scans with all of the data
combattedDf <- combattedDf[complete.cases(combattedDf), ]
print(dim(combattedDf))

# Save the resulting dataframe with combatted data and metadata
write.csv(combattedDf, fnOut, row.names = FALSE)

##
# SynthSeg ComBat

fn = '/Users/youngjm/Data/slip/fs6_stats/synthseg_2.0_phenotypes_cleaned.csv'
fnOut <- '/Users/youngjm/Data/slip/fs6_stats/06_combatted_ss_plus_metadata.csv'

analysisDf <- read.csv(fn)

prepForCombatSynthSeg(analysisDf, "/Users/youngjm/Data/slip/fs6_stats/04_toCombat_ss")

combatCommand <- paste0("python /Users/youngjm/Projects/mpr_analysis/runNeuroHarmonize.py",
                        " -p ", fnBase, "04_toCombat_ss_phenotypes.csv",
                        " -c ", fnBase, "04_toCombat_ss_covariates.csv",
                        " -o ", fnBase, "05_ss_postCombat.csv")

# Run ComBat
system(combatCommand)

# Load the combatted dataframe
combattedDf <- loadCombattedData(analysisDf, '/Users/youngjm/Data/slip/fs6_stats/05_ss_postCombat.csv')

# Keep only scans with all of the data
combattedDf <- combattedDf[complete.cases(combattedDf), ]
print(dim(combattedDf))

# Save the resulting dataframe with combatted data and metadata
write.csv(combattedDf, fnOut, row.names = FALSE)
