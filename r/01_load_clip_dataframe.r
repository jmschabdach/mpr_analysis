gc()

source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")


#-------------------------------------------------------------------------------
# Loading and Prepping Data
#
# Before loading the data, first join the IFS and FS .csv files into a single
# .csv file using the python 
#-------------------------------------------------------------------------------


##
# SLIP paper's specific process for cleaning and organizing phenotypic data
# @param fnPhenoIn A string specifying the path to the input file with all original phenotype data
# @param fnPhenoOut A string specifying the path to the output file with all cleaned phenotype data
# @param fnScannerOut A string specifying the path to the output file with the info about the scanner ids and names
# @param fnRatings A string specifying the path to the input file containing scan rating info
# @param fnDemoOut (optional) A string specifying the output path to the file containing demographic info for analyses
cleanPhenotypeDataframe <- function(fnPhenoIn, fnPhenoOut, fnScannerOut, fnRatings, fnDemoOut=""){
  
  ratingsDf <- read.csv(fnRatings)
  masterDf <- read.csv(fnPhenoIn)
  
  if ('age_in_days' %in% colnames(masterDf)){
    names(masterDf)[names(masterDf) == 'age_in_days'] <- 'age_at_scan_days'
  }
  
  if ('subj_id' %in% colnames(masterDf)){
    names(masterDf)[names(masterDf) == 'subj_id'] <- 'patient_id'
  }
  
  # Drop any rows where neurofibromatosis is in the scan_reason_categories column
  analysisDf <- masterDf[!grepl("neurofibromatosis", masterDf$scan_reason_primary), ]
  analysisDf <- analysisDf[!grepl("neurofibromatosis", analysisDf$scan_reason_categories), ]
  analysisDf <- analysisDf[(analysisDf$confirm_neurofibromatosis != TRUE), ]
  
  # If there are any missing values in the confirm_neurofibromatosis column, FALSE
  if ('1' %in% analysisDf$confirm_neurofibromatosis){
    analysisDf <- analysisDf %>%
      mutate(confirm_neurofibromatosis = case_when(
        grepl('1', confirm_neurofibromatosis, fixed=TRUE) ~ TRUE,
        grepl('0', confirm_neurofibromatosis, fixed=TRUE) ~ FALSE
      ))
  }
  
  # Drop any patient_ids from the factor levels that are no longer in the dataframe
  analysisDf$patient_id <- droplevels(as.factor(analysisDf$patient_id))
  
  # Print the number of unique patient ids in the dataframe
  print(length(levels(analysisDf$patient_id)))
  
  # Add the primary scan reason column
  analysisDf <- addPrimaryScanReasonCol(analysisDf)
  
  # Make an age in years column from the age in days column
  analysisDf$age_in_years <- analysisDf$age_at_scan_days/365.25
  
  # Some of the columns should be factors
  toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 
                'scanner_id', 'scan_reason_primary')
  analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)
  
  # Drop any 1.5T scans
  analysisDf <- analysisDf[analysisDf$MagneticFieldStrength != "1.5",]
  
  ## We only one scan per subject
  # Sort the dataframe by patient_id and scanner_id
  analysisDf <- analysisDf[ with(analysisDf, order(analysisDf$patient_id, analysisDf$scan_id)), ]
  # Drop all but the first occurrence of each patient_id
  analysisDf <- analysisDf[!duplicated(analysisDf$patient_id), ]
  # Convert patient_id to a factor
  analysisDf$patient_id <- droplevels(as.factor(analysisDf$patient_id))
  
  ## Add the average rating for each scan into the dataframe
  # If the dataframe lacks the patient id column, parse it
  if (!"patient_id" %in% colnames(analysisDf)) {
    analysisDf <- separate(analysisDf, scan_id, c("patient_id", NA, NA, NA, NA), sep="_")
  }
  # If the dataframe lacks the session id column, parse it
  if (!"sess_id" %in% colnames(analysisDf)) {
    analysisDf <- separate(analysisDf, scan_id, c(NA, "sess_id", NA, NA, NA), sep="_", remove = FALSE)
  }
  # Filter the phenotype and analysis dataframes to get only scans that have been rated
  analysisDf <- analysisDf[(analysisDf$patient_id %in% ratingsDf$subject) & (analysisDf$sess_id %in% ratingsDf$session), ] 
  ratingsDf <- ratingsDf[(ratingsDf$subject %in% analysisDf$patient_id) & (ratingsDf$session %in% analysisDf$sess_id), ] 
  
  # Add the scan grade to the phenotype dataframe
  analysisDf$average_grade <- ratingsDf$average_grade
  
  ## This section deals with column naming differences between the 
  ##   SynthSeg and FreeSurfer data
  # If the dataframe has the column CortexVol, rename other columns accordingly
  if ('CortexVol' %in% colnames(analysisDf)){
    analysisDf$GMV <- analysisDf$CortexVol
    analysisDf$WMV <- analysisDf$CerebralWhiteMatterVol
    analysisDf$sGMV <- analysisDf$SubCortGrayVol
    analysisDf$SA <- analysisDf$CorticalSurfaceArea
    analysisDf$CT <- analysisDf$MeanCorticalThickness
  # Otherwise if the Cortex column exists, rename a single column
  } else if ('Cortex' %in% colnames(analysisDf)){
    analysisDf$GMV <- analysisDf$Cortex
  }
  
  # In both cases, the TCV phenotype is a combination of the GMV and WMV phenotypes
  analysisDf$TCV <- analysisDf$GMV + analysisDf$WMV
  
  if ('BrainSegVol' %in% colnames(analysisDf)){
    analysisDf$CSF <- analysisDf$BrainSegVol - analysisDf$BrainSegVolNotVent
  } else if ('Ventricles' %in% colnames(analysisDf)) {
    analysisDf$CSF <- analysisDf$Ventricles
  }
  
  ## Drop rows and columns meeting vertain criteria
  # Drop a subject that was problematic later
  analysisDf <- analysisDf[!(grepl("sub-HM93IPIOVZ", analysisDf$patient_id)), ]
  analysisDf <- analysisDf[!(grepl("sub-HM91VO7WC1", analysisDf$patient_id)), ]
  analysisDf <- analysisDf[!(grepl("sub-HM910VH1HZ", analysisDf$patient_id)), ]
  
  # Drop columns we don't need any more (nf columns)
  toDrop <- c("X", "confirm_neurofibromatosis")
  analysisDf <- analysisDf[ , -which(names(analysisDf) %in% toDrop)]
  
  # Drop any scans with NAs
  analysisDf <- analysisDf[complete.cases(analysisDf), ]
  
  ## Rename scanner levels based on prevalence of each scanner
  analysisDf$scanner_id <- droplevels(analysisDf$scanner_id)
  metaTable <- sort(table(analysisDf$scanner_id), decreasing=TRUE)
  for (i in c(1:length(names(metaTable)))){
    newId <- paste0("Scanner ", i)
    levels(analysisDf$scanner_id)[levels(analysisDf$scanner_id) == names(metaTable)[i]] <- newId
  }
  # Save the info about the scanner names and their original ids
  scannerLookupDf <- data.frame(scanner_id = names(metaTable), 
                                figure_id=names(sort(table(analysisDf$scanner_id), decreasing=TRUE)))
  write.csv(scannerLookupDf, fnScannerOut, row.names = FALSE)
  
  # Save one phenotype file for the demographic analysis
  if ("CT" %in% colnames(analysisDf)){
    write.csv(analysisDf, fnDemoOut, row.names = FALSE)
  }
  
  # Drop any scans with ratings less than 1
  analysisDf <- analysisDf[analysisDf$average_grade >= 1, ]
  
  write.csv(analysisDf, fnPhenoOut, row.names = FALSE)
}


#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------

baseDir <- '/Users/youngjm/Data/slip/fs6_stats/'
fnRatings <- '/Users/youngjm/Data/slip/images/qc/mpr_fs_6.0.0/aggregate_ratings.csv'
fnDemoOut <- paste0(baseDir, 'original_phenotypes_demographics.csv')

## Clean and prepare the SynthSeg data -----------------------------------------
# Specify the filenames
fnPhenoIn <- paste0(baseDir,'synthseg_2.0_phenotypes.csv')
fnPhenoOut <- paste0(baseDir, 'synthseg_2.0_phenotypes_cleaned.csv')
fnScannerOut <- paste0(baseDir, 'ss_scanner_id_lookup_table.csv')

cleanPhenotypeDataframe(fnPhenoIn, fnPhenoOut, fnScannerOut, fnRatings)


## Clean and prepare the FreeSurfer 6.0 data -----------------------------------
# Specify the filenames
fnPhenoIn <- paste0(baseDir,'fs6_reconall_structural_stats_v1.csv')
fnPhenoOut <- paste0(baseDir, 'original_phenotypes_singleScanPerSubject.csv')
fnScannerOut <- paste0(baseDir, 'fs_scanner_id_lookup_table.csv')

cleanPhenotypeDataframe(fnPhenoIn, fnPhenoOut, fnScannerOut, fnRatings, fnDemoOut)

