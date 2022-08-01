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

fn = '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-07-29_highres_nocontrast_singlescanpersubject.csv'
analysisDf <- read.csv(fn)

## Step 3: ComBat the data/Load ComBat data ------------------------------------
prepForCombat <- function(df, fnBase){
  # Move metadata to the front of the dataframe
  cols <- colnames(df)
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",                
                "MagneticFieldStrength", "scanner_id", "confirm_neurofibromatosis", 
                "rawdata_image_grade", 'fs_version', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories", 'SurfaceHoles')
  phenoCols <- setdiff(cols, metaCols)
  globalPhenoCols <- c('TotalGrayVol', 'CerebralWhiteMatterVol', 'SubCortGrayVol',
                       'eTIV', 'VentricleVolume', 'CorticalSurfaceArea', 
                       'MeanCorticalThickness', 'TCV')
  nonGlobalCols <- setdiff(phenoCols, globalPhenoCols)
  
  regionalPhenoCols <- c()
  for (c in nonGlobalCols){
    if ((grepl('lh_', c) | grepl('rh_', c)) & grepl('_grayVol', c)){
      regionalPhenoCols <- append(regionalPhenoCols, c)
    }
  }
  
  newCols <- c(metaCols, globalPhenoCols, regionalPhenoCols)
  df <- df[, newCols]
  
  # Drop any scans with NAs
  df <- df[complete.cases(df), ]
  
  # Identify the scan_id + phenotypes to combat
  toCombat <- df[, c('scan_id', globalPhenoCols, regionalPhenoCols)]
  
  # Identify covariates to harmonize on/protect
  covars <- data.frame(SITE = as.character(df$scanner_id), # <-- harmonize on the first column
                       log_age_in_years = log(df$age_in_years), # <-- protect all following columns
                       sex = df$sex,
                       reason = df$top_scan_reason_factors)
  
  # Save the data and covars dataframes to csvs
  write.csv(toCombat,paste0(fnBase, "_phenotypes_toCombat.csv"), row.names = FALSE)
  write.csv(covars,paste0(fnBase, "_covariates_toCombat.csv"), row.names = FALSE)
}

prepForCombat(analysisDf, "/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats")
# prepForCombat(highQDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc1-2")
# prepForCombat(superHighQDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc2")

# STOP HERE AND RUN COMBAT VIA PYTHON

loadCombattedData <- function(df, fn){
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",                
                "MagneticFieldStrength", "scanner_id", "confirm_neurofibromatosis", 
                "rawdata_image_grade", 'fs_version', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories", 'SurfaceHoles')
  combattedDf <- read.csv(fn)
  print(dim(df))
  print(dim(combattedDf))
  # Drop rows from the dataframe if they are not in the combattedDf
  df <- df[ df$scan_id %in% combattedDf$scan_id, ]
  combattedDf <- combattedDf[ combattedDf$scan_id %in% df$scan_id, ]
  # Drop scan_ids column from the combattedDf (duplicate)
  # combattedDf <- combattedDf[ , -which(names(combattedDf) %in% c("scan_id"))]
  combattedDf <- merge(combattedDf, df[, metaCols], by='scan_id')
  print(dim(combattedDf))
  # Add the metadata back in
  # metadataDf <- df[, metaCols]
  # combattedDf <- cbind(metadataDf, combattedDf)
  return(combattedDf)
}

# Load the combatted dataframe
combattedDf <- loadCombattedData(analysisDf, '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_combatted_covariates_removed.csv')

# Keep only scans with all of the data
combattedDf <- combattedDf[complete.cases(combattedDf), ]
print(dim(combattedDf))

# Save the resulting dataframe with combatted data and metadata
write.csv(combattedDf, '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_combatted_covariates_removed_plus_metadata.csv')
