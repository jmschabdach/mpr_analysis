
# fn <- '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-07-29_highres_nocontrast_singlescanpersubject.csv'
# outFn <- "/Users/youngjm/Data/clip/fs6_stats/precentiles_rawdata_plus_metadata.csv"
# fn <- '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_combatted_plus_metadata.csv'
fn <- '/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
# outFn <- "/Users/youngjm/Data/clip/fs6_stats/precentiles_combatted_plus_metadata.csv"
outFn <- "/Users/youngjm/Data/clip/fs6_stats/precentiles_original.csv"
analysisDf <- read.csv(fn)

# Convert some columns to factors
toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 'scanner_id', 'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

## Step 4: Prep data for centilization ---------

prepForCentilizing <- function(df, fn){
  toCentilizeCols <- c('scan_id', 'age_in_years', 'sex', 'TotalGrayVol', 
                       'CerebralWhiteMatterVol', 'SubCortGrayVol', 'VentricleVolume',
                       'CorticalSurfaceArea', 'MeanCorticalThickness', 'TCV',
                       'SurfaceHoles', 'top_scan_reason_factors', 'fs_version')
  newNames <- c('participant', 'Age', 'sex', 'GMV', 'WMV', 'sGMV', 'Ventricles', 
                'SA', 'CT', 'TCV',
                'SurfaceHoles', 'top_scan_reason_factors', 'fs_version')
  toCentilizeDf <- df[, toCentilizeCols]
  toCentilizeDf <- setnames(toCentilizeDf, toCentilizeCols, newNames)
  
  toCentilizeDf <- toCentilizeDf %>%
    mutate(sex = case_when(
      grepl('F', sex, fixed=TRUE) ~ "Female",
      grepl('M', sex, fixed=TRUE) ~ "Male"
    ))
  
  toCentilizeDf$age_days <- df$age_at_scan_days + 280 # add gestational age
  toCentilizeDf$study <- paste0("JMY_CLIP_", df$scanner_id)
  print(toCentilizeDf$study)
  # toCentilizeDf$fs_version <- "FS6_T1"
  toCentilizeDf <- toCentilizeDf %>%
    mutate(fs_version = case_when(
      grepl('IFS_6.0.0', fs_version, fixed=TRUE) ~ "FSInfant",
      grepl('FS_6.0.0', fs_version, fixed=TRUE) ~ "FS6_T1"
    ))
  toCentilizeDf$country <- "USA"
  toCentilizeDf$run <- 1
  toCentilizeDf$session <- 1
  toCentilizeDf$dx <- "CN"
  reorderCols <- c('participant', 'Age', 'age_days', 'sex', 'study', 'fs_version',
                   'country', 'run', 'session', 'dx', 'top_scan_reason_factors',
                   'SurfaceHoles',
                   'GMV', 'WMV', 'sGMV', 'Ventricles', 'SA', 'CT', 'TCV')
  toCentilizeDf <- toCentilizeDf[, reorderCols]
  toCentilizeDf$INDEX.TYPE <- "NA"
  toCentilizeDf$INDEX.OB <- "NA"
  
  write.csv(toCentilizeDf, fn)
  return(toCentilizeDf)
}

preCenDf <- prepForCentilizing(analysisDf, outFn)
