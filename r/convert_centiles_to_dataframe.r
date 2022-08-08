gc()

library(plyr)
library(dplyr)

## Step 5: Load centilized phenotypes ---------------------------------

loadCentilizedFeature <- function(centileFn, phenotype){
  centileDf <- read.csv(centileFn)
  # print(head(centileDf))
  
  c1 <- paste0(phenotype, "Transformed.normalised")
  c2 <- paste0(phenotype, "Transformed.q.wre")
  # print(c(c1, c2))
  idx1 <- grep(c1, colnames(centileDf))[1]
  idx2 <- grep(c2, colnames(centileDf))[1]
  cols <- c("participant", colnames(centileDf)[idx1], colnames(centileDf[idx2]))
  
  # Get the columns of interest
  tmpDf <- centileDf[cols]
  
  # Rename two columns for joining
  colnames(tmpDf)[1] <- "scan_id"

  return(tmpDf)
}

loadCentilizedData <- function(centileDf, fnsBase){
  centileCols <- c('GMV', 'WMV', 'sGMV', 
                   'Ventricles', 'totalSA2', 'meanCT2',
                   'TCV')
  
  centileFns <- c(paste0(fnsBase, '_GMV_centiles.csv'),
                  paste0(fnsBase, '_WMV_centiles.csv'),
                  paste0(fnsBase, '_sGMV_centiles.csv'),
                  paste0(fnsBase, '_Ventricles_centiles.csv'),
                  paste0(fnsBase, '_SA_centiles.csv'),
                  paste0(fnsBase, '_CorticalThickness_centiles.csv'),
                  paste0(fnsBase, '_TotalCerebralVolume_centiles.csv'))
  
  for (i in (1:length(centileFns))){
    # centileDf <- cbind(centileDf, loadCentilizedFeature(centileFns[[i]], paste0(centileCols[[i]], 'Transformed.q.wre')))
    tmpDf <- loadCentilizedFeature(centileFns[[i]], centileCols[[i]])

    # Rename patient id to scan id
    # names(tmpDf)[names(tmpDf) == 'patient_id'] <- 'scan_id'

    # If the number of rows in the data frames is not consistent, drop rows not in both
    if (dim(centileDf)[1] != dim(tmpDf)[1]){
      centileDf <- centileDf[(centileDf$scan_id %in% tmpDf$scan_id),]
      tmpDf <- tmpDf[(tmpDf$scan_id %in% centileDf$scan_id),]
    }
    
    # # Order the data by the name of the scan
    # tmpDf <- tmpDf[order(tmpDf$scan_id), ]
    # centileDf <- centileDf[order(centileDf$scan_id), ]
    
    # centileDf <- join(centileDf, tmpDf)
    centileDf <- merge(centileDf, tmpDf, by='scan_id', all=TRUE, SORT=TRUE)
    # print(head(centileDf[colnames(tmpDf)]))
    print(dim(centileDf)) 
    
    names(centileDf)[names(centileDf) == 'centile'] <- paste0(centileCols[[i]], '_centile')
  }
  
  # Reorder scan reason factor levels
  centileDf$top_scan_reason_factors <- as.factor(centileDf$top_scan_reason_factors)
  centileDf$top_scan_reason_factors <- relevel(centileDf$top_scan_reason_factors, "other")
  
  return(centileDf)
}

#-------------------------------------------------------------------------------
# Do Stuff
#-------------------------------------------------------------------------------

# Variables to set
# fn <- '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_combatted_covariates_preserved_plus_metadata.csv'
# centileBase <- '/Users/youngjm/Data/clip/fs6_stats/combatted_covariate_effects_preserved'
# outFn <- '/Users/youngjm/Data/clip/fs6_stats/combatted_covariates_preserved_metadata_centiles.csv'

fn <- '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats.csv'
centileBase <- '/Users/youngjm/Data/clip/fs6_stats/rawdata'
outFn <- '/Users/youngjm/Data/clip/fs6_stats/combatted_rawdata_metadata_centiles.csv'

# Load the original data
analysisDf <- read.csv(fn)

if (! "top_scan_reason_factors" %in% colnames(analysisDf)){
  source("lib_mpr_analysis.r")
  analysisDf <- addPrimaryScanReasonCol(analysisDf)
}

# Convert some columns to factors
toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 'scanner_id', 'scan_reason_primary', 'top_scan_reason_factors')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

# Load the centiles
ccDf <- loadCentilizedData(analysisDf, centileBase)

# Save the new dataframe with the original phenotypes and their centiles
write.csv(ccDf, outFn)
