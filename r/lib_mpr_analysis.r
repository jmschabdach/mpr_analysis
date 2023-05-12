#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
#-------------------------------------------------------------------------------

##
#
#
#
predictCentilesForAgeRange <- function(gamModel, ageRange, euler=0, cent=0.5){
  if (euler == 0){
    newDataM <- data.frame(logAge=sort(ageRange),
                           sex=c(rep(as.factor("M"),  length(ageRange))))
    newDataF <- data.frame(logAge=sort(ageRange),
                           sex=c(rep(as.factor("F"),  length(ageRange))))
  } else  {
    print("euler!")
    newDataM <- data.frame(logAge=sort(ageRange),
                           SurfaceHoles=c(rep(euler, length(ageRange))),
                           sex=c(rep(as.factor("M"),  length(ageRange))))
    
    newDataF <- data.frame(logAge=sort(ageRange),
                           SurfaceHoles=c(rep(euler, length(ageRange))),
                           sex=c(rep(as.factor("F"),  length(ageRange))))
  } 
  
  print(cent)
  
  # Predict phenotype values for set age range for each sex
  gammModelM <- predictAll(gamModel, newdata=newDataM)
  gammModelF <- predictAll(gamModel, newdata=newDataF)
  
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
  phenoMedianPreds <- (phenoMedianPredsF + phenoMedianPredsM)/2
  
  # Return the calculated centile curve
  return(phenoMedianPreds)
}


##
# Add a column to the dataframe - CONFIRMED
addPrimaryScanReasonCol <- function(df){
  # if (!'top_scan_reason_factors' %in% colnames(df)){
  # get the top 6 primary reasons
  topScanReasons <- names(sort(table(df$scan_reason_primary), decreasing=TRUE)[1:4])
  print(topScanReasons)
    
  # build a column with these top 5 primary reasons and a generous other category
  df <- mutate(df, top_scan_reason_factors = if_else(is.element(scan_reason_primary, topScanReasons),
                                                     paste(scan_reason_primary),
                                                     "other"))
  # }
  
  df$top_scan_reason_factors <- as.factor(df$top_scan_reason_factors)
  # Put the "other" category first
  df$top_scan_reason_factors <- relevel(df$top_scan_reason_factors, "other")
  
  return(df)
}

# PULLED FROM SCRIPTS
# Calculating the centile for a subject
calculatePhenotypeCentile <- function(model, measuredPhenotypeValue, logAge, sex, surfaceHoles=c()){
  centileDistribution <- 1:9999/10000
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
    predModel <- predictAll(model, newdata=newData)
    expectedPhenotypeValue <- qGG(centileDistribution, mu=predModel$mu, sigma=predModel$sigma, nu=predModel$nu)
    centiles[i] <- centileDistribution[which.min(abs(measuredPhenotypeValue[[i]] - expectedPhenotypeValue))]
  }
  return(centiles)
}

convertAgeToYears <- function(ages){
  for (i in names(ages)){
    print(i)
    print((ages[[i]]-280)/365.25)
  }
}

calculatePipelineQc <- function(metrics1, metrics2){
  qc <- c()
  for (i in c(1:length(metrics1))) {
    qc <- append(qc, abs(metrics1[i]-metrics2[i])/((metrics1[i]+metrics2[i])/2))
  }
  return(qc)
}

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


