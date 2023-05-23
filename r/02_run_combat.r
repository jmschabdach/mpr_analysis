gc()

setwd("/Users/youngjm/Projects/mpr_analysis/")
source("r/lib_mpr_analysis.r")

#-------------------------------------------------------------------------------
# Variables to set
#-------------------------------------------------------------------------------
 
# Filepaths
baseDir <- '/Users/youngjm/Data/slip/fs6_stats/'

fnFreeSurfer <- paste0(baseDir, 'original_phenotypes_singleScanPerSubject.csv')
fnCombattedFreeSurfer <- paste0(baseDir, '06_combatted_fs_plus_metadata.csv')

fnSynthSeg <- paste0(baseDir,'synthseg_2.0_phenotypes_cleaned.csv')
fnCombattedSynthSeg <- paste0(baseDir, '06_combatted_ss_plus_metadata.csv')



#-------------------------------------------------------------------------------
# FreeSurfer ComBat
#-------------------------------------------------------------------------------

# Load the FreeSurfer dataframe
dfFreeSurfer <- read.csv(fnFreeSurfer)

# Save a subset of the FreeSurfer dataframe to a new file with ordered columns
prepForCombat(dfFreeSurfer, paste0(baseDir, "04_toCombat_fs"))

# String representing the python command to run ComBat on the prepared FreeSurfer data
combatCommand <- paste0("python ", getwd(), "/runNeuroHarmonize.py",
                        " -p ", baseDir, "04_toCombat_fs_phenotypes.csv",
                        " -c ", baseDir, "04_toCombat_fs_covariates.csv",
                        " -o ", baseDir, "05_fs_postCombat.csv")

# Run ComBat
system(combatCommand)

# Load the ComBatted FreeSurfer dataframe
dfCombattedFreeSurfer <- loadCombattedData(dfFreeSurfer, paste0(baseDir, '05_fs_postCombat.csv'))

# Keep only scans with all of the data
dfCombattedFreeSurfer <- dfCombattedFreeSurfer[complete.cases(dfCombattedFreeSurfer), ]

# Save the ComBatted FreeSurfer phenotypes and the metadata
write.csv(dfCombattedFreeSurfer, fnCombattedFreeSurfer, row.names = FALSE)

#-------------------------------------------------------------------------------
# SynthSeg ComBat
#-------------------------------------------------------------------------------

# Load the SynthSeg dataframe
dfSynthSeg <- read.csv(fnSynthSeg)

# Save a subset of the SynthSeg dataframe to a new file with ordered columns
prepForCombatSynthSeg(dfSynthSeg, paste0(baseDir, "04_toCombat_ss"))

# String representing the python command to run ComBat on the prepared SynthSeg data
combatCommand <- paste0("python ", getwd(), "/runNeuroHarmonize.py",
                        " -p ", baseDir, "04_toCombat_ss_phenotypes.csv",
                        " -c ", baseDir, "04_toCombat_ss_covariates.csv",
                        " -o ", baseDir, "05_ss_postCombat.csv")

# Run ComBat
system(combatCommand)

# Load the ComBatted SynthSeg dataframe
dfCombattedSynthSeg <- loadCombattedData(dfSynthSeg, paste0(baseDir, '05_ss_postCombat.csv'))

# Keep only scans with all of the data
dfCombattedSynthSeg <- dfCombattedSynthSeg[complete.cases(dfCombattedSynthSeg), ]
print(dim(dfCombattedSynthSeg))

# Save the resulting dataframe with combatted data and metadata
write.csv(dfCombattedSynthSeg, fnCombattedSynthSeg, row.names = FALSE)
