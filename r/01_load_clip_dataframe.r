gc()

setwd("/Users/youngjm/Projects/mpr_analysis/r/")
source("lib_mpr_analysis.r")


#-------------------------------------------------------------------------------
# Loading and Prepping Data
#
# Before loading the data, first join the IFS and FS .csv files into a single
# .csv file using the python 
#-------------------------------------------------------------------------------

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

