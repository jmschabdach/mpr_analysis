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

source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")


# globalPhenotypes <- c('TotalBrainVol', 'TotalGrayVol', 'CerebralWhiteMatterVol', 
#                       'VentricleVolume', 'SubCortGrayVol', 'CorticalSurfaceArea',
#                       'MeanCorticalThickness')

#-------------------------------------------------------------------------------
# Loading and Prepping Data
#
# Before loading the data, first join the IFS and FS .csv files into a single
# .csv file using the python 
#-------------------------------------------------------------------------------

## Step 1: load data and prep it -----------------------------------------------

# Load the dataframe containing subject demographics, imaging phenotypes, etc.
inFn <- '/Users/youngjm/Data/clip/fs6_stats/fs6_reconall_structural_stats_v1.csv'
fnOut <- '/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
scannerFnOut <- '/Users/youngjm/Data/clip/fs6_stats/fs_scanner_id_lookup_table.csv'
demoOut <- '/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_demographics.csv'

# inFn <- '/Users/youngjm/Data/clip/fs6_stats/synthseg_2.0_phenotypes.csv'
# fnOut <- '/Users/youngjm/Data/clip/fs6_stats/synthseg_2.0_phenotypes_cleaned.csv'
# scannerFnOut <- '/Users/youngjm/Data/clip/fs6_stats/ss_scanner_id_lookup_table.csv'


ratingsFn <- '/Users/youngjm/Data/clip/images/qc/mpr_fs_6.0.0/aggregate_ratings.csv'
ratingsDf <- read.csv(ratingsFn)

masterDf <- read.csv(inFn)

if ('age_in_days' %in% colnames(masterDf)){
  names(masterDf)[names(masterDf) == 'age_in_days'] <- 'age_at_scan_days'
}

if ('subj_id' %in% colnames(masterDf)){
  names(masterDf)[names(masterDf) == 'subj_id'] <- 'patient_id'
}


# Drop any rows where neurofibromatosis is in the scan_reason_categories column
analysisDf <- masterDf[!grepl("neurofibromatosis", masterDf$scan_reason_primary), ]
analysisDf <- analysisDf[!grepl("neurofibromatosis", analysisDf$scan_reason_categories), ]
# UNCOMMENT - THIS WAS REMOVED FOR SYNTHSEG ONLY
# Need to convert missing values in "confirm_neurofibromatosis" to FALSE
analysisDf <- analysisDf %>%
  mutate(confirm_neurofibromatosis = case_when(
    grepl('1', confirm_neurofibromatosis, fixed=TRUE) ~ TRUE,
    grepl('0', confirm_neurofibromatosis, fixed=TRUE) ~ FALSE
))
analysisDf <- analysisDf[(analysisDf$confirm_neurofibromatosis != TRUE), ]

# Add the primary scan reason column
analysisDf <- addPrimaryScanReasonCol(analysisDf)

# Make an age in years column from the age in days column
analysisDf$age_in_years <- analysisDf$age_at_scan_days/365.25

# Some of the columns in this df should be factors
toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 'scanner_id', 'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

# Drop any 1.5T scans
analysisDf <- analysisDf[analysisDf$MagneticFieldStrength != "1.5",]

# We only one scan per subject
# Sort the dataframe by patient_id and scanner_id
analysisDf <- analysisDf[ with(analysisDf, order(analysisDf$patient_id, analysisDf$scan_id)), ]
# Drop all but the first occurrence of each patient_id
analysisDf <- analysisDf[!duplicated(analysisDf$patient_id), ]
# Convert patient_id to a factor - idk why I did this?
analysisDf$patient_id <- droplevels(as.factor(analysisDf$patient_id))

# Incorporate the ratingsDf and the analysisDf
if (!"patient_id" %in% colnames(analysisDf)) {
  analysisDf <- separate(analysisDf, scan_id, c("patient_id", NA, NA, NA, NA), sep="_")
}
if (!"sess_id" %in% colnames(analysisDf)) {
  analysisDf <- separate(analysisDf, scan_id, c(NA, "sess_id", NA, NA, NA), sep="_", remove = FALSE)
}
analysisDf <- analysisDf[(analysisDf$patient_id %in% ratingsDf$subject) & (analysisDf$sess_id %in% ratingsDf$session), ] 
ratingsDf <- ratingsDf[(ratingsDf$subject %in% analysisDf$patient_id) & (ratingsDf$session %in% analysisDf$sess_id), ] 

analysisDf$average_grade <- ratingsDf$average_grade

# ggplot() +
#   geom_histogram(aes(x=analysisDf$average_grade), binwidth = 0.05)

# Add a column for TCV (Total Cerebrum Volume)
if (!'TCV' %in% colnames(analysisDf)){
  analysisDf$TCV <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol
}
# Add a column: TotalBrainVolume = TotalGrayVol + CerebralWhiteMatterVol + VentricleVolume + SubCortGrayVol
if (!'TCV' %in% colnames(analysisDf)) {
  if ('TotalBrainVol' %in% colnames(analysisDf)){
    names(analysisDf)[names(analysisDf) == 'TotalBrainVol'] <- 'TCV'
  } else {
    analysisDf$TCV <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol + analysisDf$VentricleVolume + analysisDf$SubCortGrayVol
  } 
}

# Drop a subject that was problematic later
analysisDf <- analysisDf[!(grepl("sub-HM93IPIOVZ", analysisDf$patient_id)), ]

# Drop columns we don't need any more (nf columns)
toDrop <- c("X", "confirm_neurofibromatosis")
analysisDf <- analysisDf[ , -which(names(analysisDf) %in% toDrop)]

# Drop any scans with NAs
analysisDf <- analysisDf[complete.cases(analysisDf), ]


# Rename scanner levels based on prevalence of each scanner
analysisDf$scanner_id <- droplevels(analysisDf$scanner_id)
print(levels(analysisDf$scanner_id))
metaTable <- sort(table(analysisDf$scanner_id), decreasing=TRUE)
for (i in c(1:length(names(metaTable)))){
  newId <- paste0("Scanner ", i)
  levels(analysisDf$scanner_id)[levels(analysisDf$scanner_id) == names(metaTable)[i]] <- newId
}
print(levels(analysisDf$scanner_id))
scannerLookupDf <- data.frame(scanner_id = names(metaTable), 
                              figure_id=names(sort(table(analysisDf$scanner_id), decreasing=TRUE)))
write.csv(scannerLookupDf, scannerFnOut, row.names = FALSE)

# Drop any scans with ratings less than 1
write.csv(analysisDf, demoOut, row.names = FALSE)
analysisDf <- analysisDf[analysisDf$average_grade >= 1, ]

write.csv(analysisDf, fnOut, row.names = FALSE)

