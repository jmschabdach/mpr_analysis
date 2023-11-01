library(ggplot2)
library(tidyr) # for drop_na
library(dplyr) # for joins
library(data.table) # for setorder
library(patchwork) # graph organization within a figure
library(gamlss) #to fit model
library(mgcv) # helps with the gam models
library(tidymv) # helps with the gam models
library(rstatix) # t-test
library(ggpubr) # show the t-test on the plot

#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS: needed for build_your_own_growth_chart.r
#-------------------------------------------------------------------------------

##
# Estimate the Nth centile for a GAMLSS model given an age range
# @param gamModel A gamlss model
# @param ageRange A list of log(post conception age in days) 
# @param euler The median SurfaceHoles value (option, produced by FreeSurfer)
# @param cent The desired centile in decimal form (default = 0.5)
predictCentilesForAgeRange <- function(gamModel, ageRange, euler=0, cent=0.5){
  if (euler == 0){
    newDataM <- data.frame(logAge=sort(ageRange),
                           sex=c(rep(as.factor("M"),  length(ageRange))))
    newDataF <- data.frame(logAge=sort(ageRange),
                           sex=c(rep(as.factor("F"),  length(ageRange))))
  } else  {
    newDataM <- data.frame(logAge=sort(ageRange),
                           SurfaceHoles=c(rep(euler, length(ageRange))),
                           sex=c(rep(as.factor("M"),  length(ageRange))))
    
    newDataF <- data.frame(logAge=sort(ageRange),
                           SurfaceHoles=c(rep(euler, length(ageRange))),
                           sex=c(rep(as.factor("F"),  length(ageRange))))
  } 
  
  # Predict phenotype values for set age range for each sex
  gammModelM <- predictAll(gamModel, newdata=newDataM, type="response")
  gammModelF <- predictAll(gamModel, newdata=newDataF, type="response")
  
  # Calculate the `cent`th centiles for the sex models
  phenoMedianPredsM <- qGG(c(cent), 
                           mu=gammModelM$mu, 
                           sigma=gammModelM$sigma, 
                           nu=gammModelM$nu)
  
  phenoMedianPredsF <- qGG(c(cent), 
                           mu=gammModelF$mu, 
                           sigma=gammModelF$sigma, 
                           nu=gammModelF$nu)
  # Average the two calculated centile curves 
  # For visualizations and correlations with the LBCC data
  phenoMedianPreds <- (phenoMedianPredsF + phenoMedianPredsM)/2
  
  # Return the calculated centile curve
  return(phenoMedianPreds)
}

##
# Estimate the centile closest to each phenotype value in a list
# @param model A GAMLSS model
# @param measuredPhenotypeValue A list of phenotype values
# @param logAge A list of log(post conception age in days, base=10)
# @param sex A list of the sex of each subject where each value is a factor
# @param surfaceHoles A list of SurfaceHoles (optional, produced by FreeSurfer)
calculatePhenotypeCentile <- function(model, measuredPhenotypeValue, logAge, sex, surfaceHoles=c()){
  centiles <- c()
  for (i in 1:length(measuredPhenotypeValue)){
    if (length(surfaceHoles) == 0) {
      newData <- data.frame(logAge=logAge[[i]],
                            sex=sex[[i]])
    } else {
      newData <- data.frame(logAge=logAge[[i]],
                            SurfaceHoles=surfaceHoles[[i]],
                            sex=sex[[i]])
    }
    predModel <- predictAll(model, newdata=newData, type="response")
    centiles[i] <- pGG(measuredPhenotypeValue[[i]], mu=predModel$mu, sigma=predModel$sigma, nu=predModel$nu)
  }
  return(centiles)
}

#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS: needed for SLIP analysis
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

##
# Add a column to the dataframe for the top reasons for a scan using a previously existing scan reason column
# @param df A dataframe containing a column called "scan_reason_primary"
addPrimaryScanReasonCol <- function(df){
  # get the top 4 primary reasons
  topScanReasons <- names(sort(table(df$scan_reason_primary), decreasing=TRUE)[1:4])
  print(topScanReasons)
    
  # build a column with these top 5 primary reasons and a generous other category
  df <- mutate(df, top_scan_reason_factors = if_else(is.element(scan_reason_primary, topScanReasons),
                                                     paste(scan_reason_primary),
                                                     "other"))
  
  df$top_scan_reason_factors <- as.factor(df$top_scan_reason_factors)
  # Put the "other" category first
  df$top_scan_reason_factors <- relevel(df$top_scan_reason_factors, "other")
  
  return(df)
}

##
# Convert a list of post conception ages from days to years
# @param ages A list of post conception ages in days
convertAgeToYears <- function(ages){
  for (i in names(ages)){
    print(i)
    print((ages[[i]]-280)/365.25)
  }
}

##
# Calculate a QC metric
# @param metrics1 A list of measures from one source
# @param metrics2 A list of the same type of measures from a different source
calculatePipelineQc <- function(metrics1, metrics2){
  qc <- c()
  for (i in c(1:length(metrics1))) {
    qc <- append(qc, abs(metrics1[i]-metrics2[i])/((metrics1[i]+metrics2[i])/2))
  }
  return(qc)
}

##
# Prepare FreeSurfer data to undergo ComBat via neuroHarmonize.py
# @param df A dataframe
# @param fnBase The path to the directory where the prepared files will be saved
prepForCombat <- function(df, fnBase){
  # Move metadata to the front of the dataframe
  cols <- colnames(df)
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",         
                "proc_ord_year",
                "MagneticFieldStrength", "scanner_id", #"confirm_neurofibromatosis", 
                "rawdata_image_grade", 'average_grade', 'fs_version', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories", 'SurfaceHoles')
  phenoCols <- setdiff(cols, metaCols)
  globalPhenoCols <- c('GMV', 'WMV', 'sGMV', 'CSF', 'TCV')
  nonGlobalCols <- setdiff(phenoCols, globalPhenoCols)
  
  regionalPhenoCols <- c()
  for (c in nonGlobalCols){
    if ((grepl('lh_', c) | grepl('rh_', c)) & grepl('_grayVol', c)){
      regionalPhenoCols <- append(regionalPhenoCols, c)
    }
  }
  
  newCols <- c(metaCols, globalPhenoCols, regionalPhenoCols)
  # print(setdiff(newCols, colnames(df)))
  df <- df[, newCols]
  
  # Drop any scans with NAs
  df <- df[complete.cases(df), ]
  print(colnames(df))
  
  # Identify the scan_id + phenotypes to combat
  toCombat <- df[, c('scan_id', globalPhenoCols, regionalPhenoCols)]
  
  # Identify covariates to harmonize on/protect
  covars <- data.frame(SITE = as.character(df$scanner_id), # <-- harmonize on the first column
                       log_age_in_years = log(df$age_in_years), # <-- protect all following columns
                       sex = df$sex,
                       SurfaceHoles = df$SurfaceHoles,
                       reason = df$top_scan_reason_factors)
  
  # Save the data and covars dataframes to csvs
  write.csv(toCombat,paste0(fnBase, "_phenotypes.csv"), row.names = FALSE)
  write.csv(covars,paste0(fnBase, "_covariates.csv"), row.names = FALSE)
}

##
# Prepare SynthSeg data to undergo ComBat via neuroHarmonize.py
# @param df A dataframe
# @param fnBase The path to the directory where the prepared files will be saved
prepForCombatSynthSeg <- function(df, fnBase){
  # Move metadata to the front of the dataframe
  cols <- colnames(df)
  print(cols)
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",     
                "proc_ord_year",
                "MagneticFieldStrength", "scanner_id", #"confirm_neurofibromatosis", 
                "rawdata_image_grade", 'average_grade', 'fs_version', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories")
  phenoCols <- setdiff(cols, metaCols)
  globalPhenoCols <- c('GMV', 'WMV', 'sGMV', 'CSF', 'TCV')
  nonGlobalCols <- setdiff(phenoCols, globalPhenoCols)
  
  regionalPhenoCols <- c()
  
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
                       # SurfaceHoles = df$SurfaceHoles,
                       reason = df$top_scan_reason_factors)
  
  # Save the data and covars dataframes to csvs
  write.csv(toCombat,paste0(fnBase, "_phenotypes.csv"), row.names = FALSE)
  write.csv(covars,paste0(fnBase, "_covariates.csv"), row.names = FALSE)
}

##
# Load the ComBatted data 
# @param df A dataframe containing demographic/metadata to be added to the ComBatted data
# @param fn The path to the file of ComBatted data to load
loadCombattedData <- function(df, fn){
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",      
                "proc_ord_year", 
                "MagneticFieldStrength", "scanner_id", #"confirm_neurofibromatosis", 
                "rawdata_image_grade", 'fs_version', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories")
  if ("SurfaceHoles" %in% colnames(df)){
    metaCols <- append(metaCols, "SurfaceHoles")
  } 
  combattedDf <- read.csv(fn)
  combattedDf$scan_id <- as.factor(combattedDf$scan_id)
  df$scan_id <- as.factor(df$scan_id)
  print(dim(df))
  print(dim(combattedDf))
  # Drop rows from the dataframe if they are not in the combattedDf
  df <- df[ df$scan_id %in% combattedDf$scan_id, ]
  print(dim(df))
  combattedDf <- combattedDf[ combattedDf$scan_id %in% df$scan_id, ]
  print(dim(combattedDf))
  # Drop scan_ids column from the combattedDf (duplicate)
  combattedDf <- merge(combattedDf, df[, metaCols], by='scan_id')
  print(dim(combattedDf))
  # Add the metadata back in
  return(combattedDf)
}


