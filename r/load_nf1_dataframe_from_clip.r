gc()

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

#-------------------------------------------------------------------------------
# Loading and Prepping Data
#
# Before loading the data, first join the IFS and FS .csv files into a single
# .csv file using the python 
#-------------------------------------------------------------------------------

# Load the dataframe containing subject demographics, imaging phenotypes, etc.
inFn <- '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats.csv'
masterDf <- read.csv(inFn)

# GRAB any rows where neurofibromatosis is in the scan_reason_categories column
analysisDf <- masterDf[grepl("neurofibromatosis", masterDf$scan_reason_primary), ]
analysisDf <- analysisDf[grepl("neurofibromatosis", analysisDf$scan_reason_categories), ]
# Need to convert missing values in "confirm_neurofibromatosis" to FALSE
analysisDf <- analysisDf %>%
  mutate(confirm_neurofibromatosis = case_when(
    grepl('1', confirm_neurofibromatosis, fixed=TRUE) ~ TRUE,
    grepl('0', confirm_neurofibromatosis, fixed=TRUE) ~ FALSE
  ))
analysisDf <- analysisDf[(analysisDf$confirm_neurofibromatosis == TRUE), ]

# Add the primary scan reason column - not relevant, scan reason is NF1
print(analysisDf$scan_reason_primary)
# analysisDf <- addPrimaryScanReasonCol(analysisDf)

# Make an age in years column from the age in days column
analysisDf$age_in_years <- analysisDf$age_at_scan_days/365.25

# Some of the columns in this df should be factors
toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 'scanner_id', 'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

# Drop any 1.5T scans
analysisDf <- analysisDf[analysisDf$MagneticFieldStrength != "1.5",]

# Drop any scans with ratings less than 0 (-1 rated scans were post contrast in Jenna's initial manual qc)
analysisDf <- analysisDf[analysisDf$rawdata_image_grade >= 0, ]

# Drop any scans with NAs
analysisDf <- analysisDf[complete.cases(analysisDf), ]

# # We only one scan per subject - not true for NF1
# # Sort the dataframe by patient_id and scanner_id
# analysisDf <- analysisDf[ with(analysisDf, order(analysisDf$patient_id, analysisDf$age_at_scan_days)), ]
# # Drop all but the first occurrence of each patient_id
# analysisDf <- analysisDf[!duplicated(analysisDf$patient_id, analysisDf$age_at_scan_days), ]

# Convert patient_id to a factor - idk why I did this?
analysisDf$patient_id <- droplevels(as.factor(analysisDf$patient_id))

# Add a column for TCV (Total Cerebrum Volume)
analysisDf$TCV <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol
# Add a column: TotalBrainVolume = TotalGrayVol + CerebralWhiteMatterVol + VentricleVolume + SubCortGrayVol
analysisDf$TotalBrainVol <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol + analysisDf$VentricleVolume + analysisDf$SubCortGrayVol

write.csv(analysisDf, '/Users/youngjm/Data/chop-nf1/fs6_stats/2022-07-29_highres_nocontrast_nf1_phenotypes.csv')
