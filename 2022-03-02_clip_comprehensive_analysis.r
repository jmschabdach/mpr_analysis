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
# library(huxtable)
# library(magrittr)
# library(MatchIt)

# Colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
#-------------------------------------------------------------------------------

generateScannerDistributionPlot <- function(dataDf, sex, colorPalette){
  # Make a set of labels for each scanner
  scannerLabels <- str_pad(c(1:length(levels(dataDf$scanner_id))), 2, pad = "0")
  scannerLabels <- paste('Scanner', scannerLabels, sep='\n')
  
  # Set the title string
  if (sex == 'M'){
    title = "Distribution of Patients Between Scanners, Male"
    dataDf <- dataDf[dataDf$sex == 'M', ]
  } else {
    title = "Distribution of Patients Between Scanners, Female"
    dataDf <- dataDf[dataDf$sex == 'F', ]
  }
  
  # dataDf <- data.frame(table(dataDf))
  
  newPlot <- ggplot(data=dataDf, aes(x=scanner_id, fill=as.factor(rawdata_image_grade))) +
    geom_bar(position="stack") +
    scale_fill_manual(values = qcColors, drop=FALSE) +
    scale_x_discrete(labels=scannerLabels, drop=FALSE) +
    labs(title = title,
         x = "Scanner ID",
         y = "# Scans per Scanner",
         fill = "Image Quality")
  
  return(newPlot)
}

generateAgeDistributionPlot <- function(dataDf, sex, colorPalette){
  # Get the max age
  xmax <- max(dataDf$age_in_years)+1
  
  # Set the title string
  if (sex == 'M'){
    title = "Distribution of Patient Age at Scan, Male"
    dataDf <- dataDf[dataDf$sex == 'M', ]
  } else {
    title = "Distribution of Patient Age at Scan, Female"
    dataDf <- dataDf[dataDf$sex == 'F', ]
  }
  newPlot <- ggplot(data=dataDf, aes(x=age_in_years, fill=as.factor(rawdata_image_grade))) +
    geom_histogram(position="stack", binwidth=91) +
    scale_fill_manual(values = qcColors, drop=FALSE) + 
    xlim(0, xmax) +
    labs(title = title,
         x = 'Patient Age (Days)',
         y = '# Patients',
         fill="Image Quality")
  
  return(newPlot)
}

##
# Given a dataframe and a column header (imaging phenotype), create a GAMM
# (GAMM = Generalized Additive Mixed Model)
# @param df A dataframe where the rows are scans and columns contain phenotypes
# @param measure Name of a column in the dataframe
# @return mixedModel The mixed model for the given formula and specificed measure
createScanReasonGamm <- function(df, measure) {
  formula <- as.formula(paste(measure, "s(log(age_in_years), fx=T) +
                      SurfaceHoles +
                      top_scan_reason_factors +
                      sex", sep="~"))
  
  # Make a basic linear model accounting for age and surface holes for the given data frame
  mixedModel <- gamm(formula,
                     random = list(scanner_id=~1),
                     data = df)
  # Return the model
  return(mixedModel)
}


##
# Make a scatter plot of a specific measure for the original data. Include the
# confidence intervals that come from the output of the GAMM predictions
# @param predictions Output from predict_gam() for the measure and its gam
# @param origData A dataframe containing the original data
# @param measure A column of the original dataframe
# @param measureTitle A string used to describe the measure in the plot
# @returns plot01 A ggplot object with scatter points and a confidence interval
generatePlotScatterWithCI <- function(predictions, origData, measure, measureTitle){
  # Colorblind palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  plot01 <- predictions %>%
    ggplot(aes(x=age_in_years, y=fit), log='x') +
    geom_point(data=origData, aes(x=age_in_years, y=measure, color=top_scan_reason_factors)) +
    geom_smooth_ci() +
    scale_color_manual(values = cbbPalette, name = "Scan Reason Category") +
    # scale_x_continuous(trans='log10') +
    theme(axis.title = element_blank()) +
    labs(title = measureTitle)
  
  print(plot01)
  
  return(plot01)
}


##
# Make a table showing the statistical significance of nonlinear parameters on image phenotype
# @param aGam A gam object
# @param nTests A number indicating how many tests are being performed/should be corrected for
# @returns paraT A table object to be displayed using gridExtra
generateScanReasonsParametricTable <- function(aGam){
  subTab <- signif(summary(aGam)$s.table, 3)
  rownames(subTab) <- c('Age')
  
  # Add the p.adjust
  # subTab <- cbind(subTab, p.adjust(subTab[,"p-value"], n = nTests))
  colnames(subTab) <- c('edf', 'Ref.df', 'F', 'p-value')
  cols <- c("p-value")
  subTab <- subTab[,cols]

  return(subTab)
}


##
# Make a table showing the statistical significance of linear parameters on image phenotype
# @param aGam A gam object
# @param nTests A number indicating how many tests are being performed/should be corrected for
# @returns linT A table object to be displayed using gridExtra
generateScanReasonsLinearTable <- function(aGam){
  subTab <- signif(summary(aGam)$p.table, 3)
  
  newRowNames <- c('Surface Holes',
                   'Scan Reason: Developmental Disorder',
                   'Scan Reason: Eye/Vision Finding',
                   'Scan Reason: Headaches',
                   'Scan Reason: Non-Brain Lesion',
                   'Scan Reason: Seizures',
                   'Sex')
  rows <- rownames(subTab)[2:length(rownames(subTab))]
  subTab <- subset(subTab, rownames(subTab) %in% rows)
  rownames(subTab) <- newRowNames

  # Rename and trim the columns
  colnames(subTab) <- c('est', 'std err', 't-value', 'p-value')
  cols <- c('p-value')
  subTab <- subTab[,cols]

  return(subTab)
}

##
# Add a column to the dataframe 
addPrimaryScanReasonCol <- function(df){
  # get the top 6 primary reasons
  topScanReasons <- names(sort(table(df$scan_reason_primary), decreasing=TRUE)[1:5])
  
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
# Make a table showing the statistical significance of linear parameters on image phenotype
# @param aGam A gam object
# @param nTests A number indicating how many tests are being performed/should be corrected for
# @returns linT A table object to be displayed using gridExtra
getBetasFromGamSummary <- function(aGam){
  subTab <- signif(summary(aGam)$p.table, 3)
  
  newRowNames <- c('Surface Holes',
                   'Scan Reason: Developmental Disorder',
                   'Scan Reason: Eye/Vision Finding',
                   'Scan Reason: Headaches',
                   'Scan Reason: Non-Brain Lesion',
                   'Scan Reason: Seizures',
                   'Sex')
  rows <- rownames(subTab)[2:length(rownames(subTab))]
  subTab <- subset(subTab, rownames(subTab) %in% rows)
  rownames(subTab) <- newRowNames
  
  # Rename and trim the columns
  colnames(subTab) <- c('Beta', 'std err', 't-value', 'p-value')
  cols <- c('Beta')
  subTab <- subTab[,cols]
  
  return(subTab)
}

getAgeAtPeak <- function(preds){
  # calculate the gam
  # Get the peak of the gam
  peak <- preds %>%
    ungroup %>%
    slice_max(fit) %>%
    pull(age_in_years)
  return(peak)
}

##
# Make the composite plot and tables for a given dataframe/set of phenotypes
# @param df A dataframe where rows are scans and columns are phenotypes and features
# @param phenotypes A vector of phenotypes to analyze
# @param colNames A vector of cleaned names to use for the table columns
# @param rowNames A vector of cleaned names to use for the table rows
# @param title A string to differentiate the figures for this dataframe vs. other dataframes
generatePlotsAndTablesForDataset <- function(df, phenotypes, colNames, rowNames, title){
  # Initialize variables
  # Set up empty variables to generate table later
  ageAtPeak <- list()
  predictionPlots <- list()
  tableValues <- list()
  tableBetas <- list()
  
  modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
                           sex='M', 
                           top_scan_reason_factors='headaches')
  
  # for phenotype in phenotypes...
  for (phenotype in phenotypes){
    # Generate GAMs using the formula that incorporates scan reason
    gammScanReason <- createScanReasonGamm(df, phenotype)
    
    # Predict on the GAMs for the actual data
    gammScanReasonPreds <- predict_gam(gammScanReason$gam, values = modelFixedValues)
    
    # Get the age at peak
    ageAtPeak[[phenotype]] <- getAgeAtPeak(gammScanReasonPreds)
    
    # Generate scatter plots with confidence intervals
    plotScatterScanReason <- generatePlotScatterWithCI(gammScanReasonPreds, df, df[ , phenotype], phenotype)
    predictionPlots[[phenotype]] <- plotScatterScanReason
    
    tParam <- generateScanReasonsParametricTable(gammScanReason$gam)
    tLin <- generateScanReasonsLinearTable(gammScanReason$gam)
    tableValues[[phenotype]] <- append(tLin, list(Age = tParam))
    tableBetas[[phenotype]] <- getBetasFromGamSummary(gammScanReason$gam)
  }
  
  # Plot all scatter/CI plots in 1 figure
  patchwork <- wrap_plots(predictionPlots, ncol=2, guides="collect")
  print(patchwork + plot_annotation(title=paste("Trajectories of", title)))
  
  # Make table: p-values of factor/phenotype
  firstStep <- lapply(tableValues, unlist) 
  secondStep <- as.data.frame(firstStep, stringsAsFactors = F) 
  colnames(secondStep) <- colNames
  rownames(secondStep) <- rowNames
  t1 <- gridExtra::tableGrob(secondStep)
  gridExtra::grid.arrange(top=paste("p-values of", title), t1)

  # Make table: betas
  firstStep <- lapply(tableBetas, unlist) 
  secondStep <- as.data.frame(firstStep, stringsAsFactors = F) 
  colnames(secondStep) <- tableColumns
  rownames(secondStep) <- tableRows[1:length(tableRows)-1]
  t1 <- gridExtra::tableGrob(secondStep)
  gridExtra::grid.arrange(top=paste("Betas of", title), t1)
  
  return(ageAtPeak)
}

#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------

##------------------------------------------------------------------------------
# Step 1: load data and prep it

# Load the master file containing subject demographics, imaging phenotypes, etc.
inFn <- '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-03_analysis_features.csv'
# inFn <- '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_7.1.1_structural_stats.csv'
masterDf <- read.csv(inFn)

# Drop any rows where neurofibromatosis is in the scan_reason_categories column
analysisDf <- masterDf[!grepl("neurofibromatosis", masterDf$scan_reason_primary), ]
analysisDf <- analysisDf[!grepl("neurofibromatosis", analysisDf$scan_reason_categories), ]
analysisDf <- analysisDf[!grepl("true", analysisDf$confirm_neurofibromatosis), ]

# Make an age in years column
analysisDf$age_in_years <- analysisDf$age_in_days/365.25

# Add new column: going to sum TotalGrayVol + CerebralWhiteMatterVol + VentricleVolume + SubCortGrayVol
analysisDf$TotalBrainVol <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol + analysisDf$VentricleVolume + analysisDf$SubCortGrayVol

# Some of the columns in this df should be factors
toFactor <- c('sex', 'Processing', 'MagneticFieldStrength', 'scanner_id', 
              'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

##------------------------------------------------------------------------------
# Step 2: generate basic demographic plots

qcColors <- c("#999999", "#0072B2", "#56B4E9")

# Make a plot for the distribution of ages with coloring for sex 
pAgesMale <- generateAgeDistributionPlot(analysisDf, 'M', qcColors)
pAgesFemale <- generateAgeDistributionPlot(analysisDf, 'F', qcColors)

pScannersMale <- generateScannerDistributionPlot(analysisDf, 'M', qcColors)
pScannersFemale <- generateScannerDistributionPlot(analysisDf, 'F', qcColors)

# Make plots for QC
pQcSurfaceHoles <- analysisDf %>%
  ggplot(aes(x=as.factor(rawdata_image_grade), y=SurfaceHoles)) +
  geom_violin(aes(fill=as.factor(rawdata_image_grade))) +
  geom_boxplot(width=0.1)+ 
  scale_fill_manual(values=qcColors) +
  labs(title = "Distribution of Surface Holes Count at Scan\nAcross Image QC Ratings",
       x = 'Image QC Rating Group',
       y = 'Number of Surface Holes (count)',
       fill="Image Quality")

pQcAge <- analysisDf %>%
  ggplot(aes(x=as.factor(rawdata_image_grade), y=age_in_years)) +
  geom_violin(aes(fill=as.factor(rawdata_image_grade))) +
  geom_boxplot(width=0.1)+ 
  scale_fill_manual(values=qcColors)+
  labs(title = "Distribution of Age at Scan\nAcross Image QC Ratings",
       x = 'Image QC Rating Group',
       y = 'Age at Scan (Days)',
       fill="Image Quality")

grob <- patchworkGrob(pAgesMale + pScannersMale + 
                        pAgesFemale + pScannersFemale +
                        pQcAge + pQcSurfaceHoles +
                        plot_layout(guides="collect", ncol = 2))
gridExtra::grid.arrange(grob)

##------------------------------------------------------------------------------
# Step 3: Evaluate QC
# [ ] 4.1 LM for age, scanner, quality?
# [x] 4.2 KS tests for scanner holes and quality

# ks.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, 
#         analysisDf[analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles)

# Compare number of surface holes (Euler number) to image grade
mylogit <- glm(as.factor(rawdata_image_grade) ~ SurfaceHoles, data=analysisDf, family='binomial')
summary(mylogit)

##------------------------------------------------------------------------------
# Drop the rows where rawdata_image_grade <= 1
analysisDf <- analysisDf[(analysisDf$rawdata_image_grade >= 1),]

##------------------------------------------------------------------------------
# Step 4: Make a plot for the distribution of reasons for a scan

analysisDf <- addPrimaryScanReasonCol(analysisDf)
reasonsTable <- summary(analysisDf$top_scan_reason_factors)
pieLabels <- paste(names(reasonsTable), reasonsTable, sep='\n')
pie(reasonsTable, pieLabels,
    main="Top Scan Reasons\n(with sample sizes)")

##------------------------------------------------------------------------------
# Step 5: Generate GAMs for each phenotype of interest

phenotypes <- c('TotalBrainVol', 'TotalGrayVol', 'CerebralWhiteMatterVol', 
                'VentricleVolume', 'SubCortGrayVol', 'SumCorticalSurfaceArea', 'AvgCorticalThickAvg')

# Build a table with all of the global phenotype p-vals parameters
tableColumns <- c("Total Brain\nVolume", "Total Gray\nVolume", "Cerebral White\nMatter Volume",
                  "Ventricle\nVolume", "Subcortical Gray\nVolume", "Total Cortical\nSurface Area",
                  "Average Cortical\nThickness")
tableRows <- c("Surface Holes", "Scan Reason:\nDevelopmental Disorder",
               "Scan Reason:\nEye/Vision Finding", "Scan Reason:\nHeadaches",
               "Scan Reason:\nNon-Brain Lesion", "Scan Reason:\nSeizures",
               "Sex = Male", "Age")

rawAgeAtPeak <- generatePlotsAndTablesForDataset(analysisDf, phenotypes, tableColumns, tableRows, 'Raw Phenotypes')


##------------------------------------------------------------------------------
# Step 6: Run the GAM analysis for the centile scores

# Load the master file containing subject demographics, imaging phenotypes, etc.
# centileFn <- '/Users/youngjm/Data/clip/tables/imaging_results/2022-03-analysis_feature_centilescores.csv'
centileFn <- '/Users/youngjm/Data/clip/tables/imaging_results/2022-03_CLIP_lifespan_centilesGMV.csv'
centileDf <- read.csv(centileFn)

# Add the age column without postgestational age in days
centileDf$age_in_years <- (centileDf$age_in_days - 280)/365.25

# Make sure that sex is M/F
centileDf$sex <- case_when(
  centileDf$sex == "Male" ~ "M",
  centileDf$sex == "Female" ~ "F"
)

# Some of the columns in this df should be factors
toFactor <- c('sex', 'study', 'fs_version')
centileDf[toFactor] <- lapply(centileDf[toFactor], factor)

# Add columns from 
centileDf$rawdata_image_grade <- masterDf$rawdata_image_grade
centileDf$scan_reason_primary <- masterDf$scan_reason_primary
centileDf$scan_reason_categories <- masterDf$scan_reason_categories
centileDf$confirm_neurofibromatosis <- masterDf$confirm_neurofibromatosis

# Drop any rows where neurofibromatosis is in the scan_reason_categories column
centileDf <- centileDf[!grepl("neurofibromatosis", centileDf$scan_reason_primary), ]
centileDf <- centileDf[!grepl("neurofibromatosis", centileDf$scan_reason_categories), ]
centileDf <- centileDf[!grepl("true", centileDf$confirm_neurofibromatosis), ]

# Remove the low quality scans
centileDf <- centileDf[(centileDf$rawdata_image_grade >= 1),]

# Get primary reason from the analysisDf
centileDf$top_scan_reason_factors <- analysisDf$top_scan_reason_factors
centileDf$SurfaceHoles <- analysisDf$SurfaceHoles
centileDf$scanner_id <- analysisDf$scanner_id

# Set up list of phenotypes to pull from the centiles
phenotypes <- c('GMV_quantile', 'WMV_quantile', 
                'Ventricles_quantile', 'sGMV_quantile') #, 'SumCorticalSurfaceArea', 'AvgCorticalThickAvg')

# Build a table with all of the global phenotype p-vals parameters
tableColumns <- c("Total Gray\nVolume (Centile)", "Cerebral White Matter\nVolume (Centile)",
                  "Ventricle Volume\n(Centile)", "Subcortical Gray\nVolume (Centile)")
tableRows <- c("Surface Holes", "Scan Reason:\nDevelopmental Disorder",
               "Scan Reason:\nEye/Vision Finding", "Scan Reason:\nHeadaches",
               "Scan Reason:\nNon-Brain Lesion", "Scan Reason:\nSeizures",
               "Sex = Male", "Age")

centileAgeAtPeak <- generatePlotsAndTablesForDataset(centileDf, phenotypes, tableColumns, tableRows, 'Centilized Phenotypes')


##------------------------------------------------------------------------------
# Step 7: ComBat the data

# # Load the ComBat library
# library(sva)
# 
# # Pull the scan ids
# scanIds <- analysisDf$scan_id
# # Pull the phenotypes we want to look at into a DF where rows are features and cols are scans
# toCombat <- as.data.frame(t(analysisDf[c('TotalBrainVol', 'TotalGrayVol', 'CerebralWhiteMatterVol', 
#                                          'VentricleVolume', 'SubCortGrayVol', 'SumCorticalSurfaceArea', 
#                                          'AvgCorticalThickAvg')]))
# # Set the column names of the isolated phenotypes
# colnames(toCombat) <- scanIds
# # We want to remove differences based on the scanner id, so pull that info to use as the batch variable
# batch <- analysisDf$scanner_id
# 
# # Run ComBat
# combattedDf <- ComBat(toCombat, batch)
# 
# # Pull the column names of the combatted data
# # scanIds <- colnames(combattedDf)
# # Transpose the combatted data
# combattedDf <- as.data.frame(t(combattedDf))
# # Add in other metadata needed for analysis
# # combattedDf$scan_id <- scanIds
# metadataCols <- colnames(analysisDf)[1:15]
# metadataDf <- analysisDf[, metadataCols]
# combattedDf <- cbind(metadataDf, combattedDf)
# combattedDf$SurfaceHoles <- analysisDf$SurfaceHoles
# combattedDf$top_scan_reason_factors <- analysisDf$top_scan_reason_factors
# combattedDf$age_in_years <- analysisDf$age_in_years
# # save the combatted data
# write.csv(combattedDf,"/Users/youngjm/Data/clip/tables/imaging_results/2022-03_combatted_global_phenotypes.csv", row.names = FALSE)

combattedDf <- read.csv("/Users/youngjm/Data/clip/tables/imaging_results/2022-03_combatted_global_phenotypes.csv")

##------------------------------------------------------------------------------
# Step 8: GAMs for the combatted data

phenotypes <- c('TotalBrainVol', 'TotalGrayVol', 'CerebralWhiteMatterVol', 
                'VentricleVolume', 'SubCortGrayVol', 'SumCorticalSurfaceArea', 'AvgCorticalThickAvg')

# Build a table with all of the global phenotype p-vals parameters
tableColumns <- c("Total Brain\nVolume", "Total Gray\nVolume", "Cerebral White\nMatter Volume",
                  "Ventricle\nVolume", "Subcortical Gray\nVolume", "Total Cortical\nSurface Area",
                  "Average Cortical\nThickness")
tableRows <- c("Surface Holes", "Scan Reason:\nDevelopmental Disorder",
               "Scan Reason:\nEye/Vision Finding", "Scan Reason:\nHeadaches",
               "Scan Reason:\nNon-Brain Lesion", "Scan Reason:\nSeizures",
               "Sex = Male", "Age")
# 1600x2400
combatAgeAtPeak <- generatePlotsAndTablesForDataset(combattedDf, phenotypes, tableColumns, tableRows, 'Phenotypes after ComBat')

##------------------------------------------------------------------------------
# Step 9: Make a table of the age at peak

ageAtPeaks <- rbind(rawAgeAtPeak, combatAgeAtPeak)
ageAtPeaks <- do.call(rbind, Map(data.frame, RawPheno=rawAgeAtPeak, 
                                 ComBatPheno=combatAgeAtPeak))

# Make table: betas
firstStep <- lapply(ageAtPeaks, unlist)
secondStep <- as.data.frame(ageAtPeaks, stringsAsFactors = F) 
rownames(secondStep) <- tableColumns
colnames(secondStep) <- c('Raw Phenotypes', 'Combatted Phenotypes')
t1 <- gridExtra::tableGrob(secondStep)
gridExtra::grid.arrange(top="Age At Peak (years)", t1)

##------------------------------------------------------------------------------
# Step 10: Regional phenotypes

metadataCols <- colnames(analysisDf)[1:15]
metadataCols <- append(metadataCols, 'age_in_years')
metadataCols <- append(metadataCols, 'SurfaceHoles')
metadataCols <- append(metadataCols, 'top_scan_reason_factors')
greyVolCols <- analysisDfCols[endsWith(analysisDfCols, '_grayVol')]
grayVolLhCols <- greyVolCols[startsWith(greyVolCols, 'lh_')]
grayVolRhCols <- greyVolCols[startsWith(greyVolCols, 'rh_')]

regionalPhenotypes <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal',
                        # 'corpus callosum', 
                        'cuneus', 'entorhinal', 'frontal pole', 
                        'fusiform', 'inferior parietal', 'inferior temporal', 'insula',
                        'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal',
                        'lingual', 'medial orbitofrontal', 'middle temporal', 'paracentral',
                        'parahippocampal', 'pars opercularis', 'pars orbitalis',
                        'pars triangularis', 'pericalcarine', 'postcentral',
                        'posterior cingulate', 'precentral', 'precuneus', 
                        'rostral anterior cingulate', 'rostral middle frontal',
                        'superior frontal', 'superior parietal', 'superior temporal',
                        'supramarginal', 'temporal pole', 'transverse temporal')

# Get the data for the volCols
grayVolLhDf <- analysisDf[, append(metadataCols, grayVolLhCols) ]
grayVolRhDf <- analysisDf[, append(metadataCols, grayVolRhCols) ]

# # Rename columns
# colnames(grayVolLhDf) <- append(metadataCols, regionalPhenotypes)
# colnames(grayVolRhDf) <- append(metadataCols, regionalPhenotypes)

##
# Make the composite plot and tables for a given dataframe/set of phenotypes
# @param df A dataframe where rows are scans and columns are phenotypes and features
# @param phenotypes A vector of phenotypes to analyze
# @param colNames A vector of cleaned names to use for the table columns
# @param rowNames A vector of cleaned names to use for the table rows
# @param title A string to differentiate the figures for this dataframe vs. other dataframes
generateHemispherePlots <- function(lhDf, rhDf, phenotypes, title){
  # Initialize variables
  # Set up empty variables to generate table later
  ageAtPeak <- c()
  hemispheres <- c()
  regions <- c()
  
  modelFixedValues <- list(SurfaceHoles = mean(lhDf$SurfaceHoles), 
                           sex='M', 
                           top_scan_reason_factors='headaches')
  
  # for phenotype in phenotypes...
  for (phenotype in phenotypes){
    print(phenotype)
    parsedPheno <- gsub(" ", "", phenotype, fixed=TRUE)
    # Generate GAMs using the formula that incorporates scan reason
    gammLh <- createScanReasonGamm(lhDf, paste('lh', parsedPheno, 'grayVol', sep='_'))
    gammRh <- createScanReasonGamm(rhDf, paste('rh', parsedPheno, 'grayVol', sep='_'))
    
    # Predict on the GAMs for the actual data
    gammLhPreds <- predict_gam(gammLh$gam, values = modelFixedValues)
    gammRhPreds <- predict_gam(gammRh$gam, values = modelFixedValues)
    
    # Get the age at peak
    ageAtPeak <- append(ageAtPeak, c(getAgeAtPeak(gammLhPreds), getAgeAtPeak(gammRhPreds)))
    hemispheres <- append(hemispheres, c('left', 'right'))
    regions <- append(regions, c(phenotype, phenotype))
  }
  
  # Convert the lists into a dataframe
  results = as.data.frame(cbind(region=regions, 
                                em=ageAtPeak, 
                                hemi=hemispheres),
                       stringsAsFactors=F)
  
  # Plot the dataframe using ggseg
  p <- results %>% 
    ggseg(mapping=aes(fill=as.numeric(em)),
        # show.legend = F,
        position ='stacked')
  
  grid.arrange(p)
}


generateHemispherePlots(grayVolLhDf, grayVolRhDf, regionalPhenotypes, '')  
