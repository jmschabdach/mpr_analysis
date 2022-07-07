gc()
dev.off(dev.list()["RStudioGD"])

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
    title = "Scanner Distribution (Male)"
    dataDf <- dataDf[dataDf$sex == 'M', ]
  } else {
    title = "Scanner Distribution (Female)"
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
         fill = "Image Quality") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(newPlot)
}

generateAgeDistributionPlot <- function(dataDf, sex, colorPalette){
  # Get the max age
  xmax <- max(dataDf$age_in_years)+1
  print(xmax)
  
  # Set the title string
  if (sex == 'M'){
    title = "Age at Scan (Male)"
    dataDf <- dataDf[dataDf$sex == 'M', ]
  } else {
    title = "Age at Scan (Female)"
    dataDf <- dataDf[dataDf$sex == 'F', ]
  }
  
  newPlot <- ggplot(data=dataDf, aes(x=age_in_years, fill=as.factor(rawdata_image_grade))) +
    geom_histogram(position=position_stack(), binwidth = 0.5) +
    scale_fill_manual(values = qcColors, drop=FALSE) + 
    xlim(0, xmax) +
    labs(title = title,
         x = 'Patient Age (Years)',
         y = '# Patients',
         fill="Image Quality") +
    theme(plot.title = element_text(hjust = 0.5))
  
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
  
  # modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
  #                          sex='M', 
  #                          top_scan_reason_factors='headaches')
  
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
  
  peakAge <- getAgeAtPeak(predictions)
  
  plot01 <- predictions %>%
    ggplot(aes(x=age_in_years, y=fit), log='x') +
    geom_point(data=origData, alpha=0.5, aes(x=age_in_years, y=measure, color=top_scan_reason_factors)) +
    geom_smooth_ci() +
    # geom_vline(xintercept=peakAge) +
    scale_color_manual(values = cbbPalette, name = "Scan Reason Category") +
    # scale_x_continuous(trans='log10') +
    theme(axis.title = element_blank()) +
    labs(title = measureTitle)
  
  print(plot01)
  
  return(plot01)
}

generateDiagnosisPlotCI <- function(predictions, origData, measure, measureTitle){
  
  dxColors <- c("#000000", "#E69F00")
  
  plot02 <- predictions %>%
    ggplot(aes(x=age_in_years, y=fit)) +
    geom_smooth_ci() +
    scale_color_manual(values = dxColors) +
    # scale_x_continuous(trans='log10') +
    theme(axis.title = element_blank()) 
    # labs(title = "Growth Trajectory Over Lifespan") #paste(measureTitle, ": Growth Trajectories"))
  
  return(plot02)
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

## Step 1: load data and prep it ------------------------------------------------

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

## Step 2: generate basic demographic plots ------------------------------------
# 

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
  labs(title = "Euler Number vs. QC Ratings",
       x = 'Image QC Rating Group',
       y = 'Number of Surface Holes (count)',
       fill="Image Quality") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

pQcAge <- analysisDf %>%
  ggplot(aes(x=as.factor(rawdata_image_grade), y=age_in_years)) +
  geom_violin(aes(fill=as.factor(rawdata_image_grade))) +
  geom_boxplot(width=0.1)+ 
  scale_fill_manual(values=qcColors)+
  labs(title = "Age at Scan vs. QC Ratings",
       x = 'Image QC Rating Group',
       y = 'Age at Scan (Days)',
       fill="Image Quality") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

grob <- patchworkGrob(pAgesMale + pAgesFemale + 
                        pScannersMale + pScannersFemale +
                        pQcAge + pQcSurfaceHoles +
                        plot_layout(guides="collect", ncol = 2))
gridExtra::grid.arrange(grob)

## Step 3: Evaluate QC ---------------------------------------------------------

# ks.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, 
#         analysisDf[analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles)

# Compare number of surface holes (Euler number) to image grade
mylogit <- glm(as.factor(rawdata_image_grade) ~ SurfaceHoles, data=analysisDf, family='binomial')
summary(mylogit)

## Drop the rows where rawdata_image_grade <= 1 --------------------------------
analysisDf <- analysisDf[(analysisDf$rawdata_image_grade >= 1),]

## Step 4: Make a plot for the distribution of reasons for a scan --------------

analysisDf <- addPrimaryScanReasonCol(analysisDf)
reasonsTable <- summary(analysisDf$top_scan_reason_factors)
pieLabels <- paste(names(reasonsTable), reasonsTable, sep='\n')
pie(reasonsTable, pieLabels,
    main="Top Scan Reasons\n(with sample sizes)")

## Step 5: Generate GAMs for each phenotype of interest ------------------------

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


## Step 6: Run the GAM analysis for the centile scores -------------------------


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


## Step 7: ComBat the data -----------------------------------------------------

# # Load the ComBat library
# library(sva)
# 
# Pull the scan ids
scanIds <- analysisDf$scan_id
# Pull the phenotypes we want to look at into a DF where rows are features and cols are scans
phenotypes <- c('TotalBrainVol', 'TotalGrayVol', 'CerebralWhiteMatterVol',
                'VentricleVolume', 'SubCortGrayVol', 'SumCorticalSurfaceArea',
                'AvgCorticalThickAvg')
# analysisDfCols <- colnames(analysisDf)
# phenotypes <- c(phenotypes, analysisDfCols[endsWith(analysisDfCols, '_grayVol')])
toCombat <- analysisDf[phenotypes]
# Set the column names of the isolated phenotypes
toCombat$scan_id <- scanIds
# We want to remove differences based on the scanner id, so pull that info to use as the batch variable
# We want to preserve differences between age, sex, dx, and FS version
covars <- data.frame(SITE = analysisDf$scanner_id, 
                    FSVersion = analysisDf$Processing,
                    age_in_years = analysisDf$age_in_years,
                    sex = analysisDf$sex,
                    reason = analysisDf$top_scan_reason_factors)

colnames(covars)
colnames(toCombat)

# Save the data and covars dataframes to csvs
write.csv(toCombat,"/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-04_clip_phenotypes_toCombat.csv", row.names = FALSE)
write.csv(covars,"/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-04_clip_covariates_toCombat.csv", row.names = FALSE)

combattedDf <- read.csv('/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-04_clip_phenotypes_combatted.csv')

# # Run ComBat
# combattedDf <- ComBat(toCombat, batch)
# 
# # Pull the column names of the combatted data
# # scanIds <- colnames(combattedDf)
# # Transpose the combatted data
# combattedDf <- as.data.frame(t(combattedDf))
# # Add in other metadata needed for analysis
# # combattedDf$scan_id <- scanIds
metadataCols <- colnames(analysisDf)[1:15]
metadataDf <- analysisDf[, metadataCols]
combattedDf <- cbind(metadataDf, combattedDf)
combattedDf$SurfaceHoles <- analysisDf$SurfaceHoles
combattedDf$top_scan_reason_factors <- analysisDf$top_scan_reason_factors
combattedDf$age_in_years <- analysisDf$age_in_years
# save the combatted data
write.csv(combattedDf,"/Users/youngjm/Data/clip/tables/imaging_results/2022-04_clip_phenotypes_combatted.csv", row.names = FALSE)

combattedDf <- read.csv("/Users/youngjm/Data/clip/tables/imaging_results/2022-04_clip_phenotypes_combatted.csv")

## Step 8: GAMs for the combatted data -----------------------------------------

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

## Step 9: Make a table of the age at peak -------------------------------------

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

## Step 10: Make ggseg figures for regional phenotypes -------------------------

# analysisDf <- analysisDf[analysisDf$rawdata_image_grade > 1,]

# brainDf <- combattedDf
brainDf <- analysisDf

metadataCols <- colnames(brainDf)[1:15]
metadataCols <- append(metadataCols, 'age_in_years')
metadataCols <- append(metadataCols, 'SurfaceHoles')
metadataCols <- append(metadataCols, 'top_scan_reason_factors')
brainDfCols <- colnames(brainDf)

regionalPhenotypes <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal',
                        # 'corpus callosum', 
                        'cuneus', # 'entorhinal', 'frontal pole', # these two were excluded from Lifespan
                        'fusiform', 'inferior parietal', 'inferior temporal', 'insula',
                        'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal',
                        'lingual', 'medial orbitofrontal', 'middle temporal', 'paracentral',
                        'parahippocampal', 'pars opercularis', 'pars orbitalis',
                        'pars triangularis', 'pericalcarine', 'postcentral',
                        'posterior cingulate', 'precentral', 'precuneus', 
                        'rostral anterior cingulate', 'rostral middle frontal',
                        'superior frontal', 'superior parietal', 'superior temporal',
                        'supramarginal', #'temporal pole', # excluded from Lifespan
                        'transverse temporal')
parsedRegionalPhenotypes <- c('bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal',
                              # 'corpus callosum', 
                              'cuneus', #'entorhinal', 'frontalpole', 
                              'fusiform', 'inferiorparietal', 'inferiortemporal', 'insula',
                              'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal',
                              'lingual', 'medialorbitofrontal', 'middletemporal', 'paracentral',
                              'parahippocampal', 'parsopercularis', 'parsorbitalis',
                              'parstriangularis', 'pericalcarine', 'postcentral',
                              'posteriorcingulate', 'precentral', 'precuneus', 
                              'rostralanteriorcingulate', 'rostralmiddlefrontal',
                              'superiorfrontal', 'superiorparietal', 'superiortemporal',
                              'supramarginal', # 'temporalpole', 
                              'transversetemporal')


greyVolCols <- brainDfCols[endsWith(brainDfCols, '_grayVol')]
print(greyVolCols)
phenoCols <- c()
for (pheno in parsedRegionalPhenotypes) {
  cols <- greyVolCols[grepl(paste("_",pheno, sep=''), greyVolCols, fixed=TRUE)]
  phenoCols <- c(phenoCols, cols)
}
grayVolLhCols <- sort(phenoCols[startsWith(phenoCols, 'lh_')])
grayVolRhCols <- sort(phenoCols[startsWith(phenoCols, 'rh_')])

# Let's calculate the bilateral average of each phenotype
for (i in 1:length(parsedRegionalPhenotypes)){
  brainDf[[parsedRegionalPhenotypes[[i]]]] <- (brainDf[[grayVolLhCols[[i]]]] + brainDf[[grayVolRhCols[[i]]]])/2
}

greyVolDf <- brainDf[append(metadataCols, parsedRegionalPhenotypes)]

# Generate GAMMS for regional phenotypes
generatePlotsAndTablesForDataset(brainDf, parsedRegionalPhenotypes, tableColumns, tableRows, 'Regional Phenotypes')

##
# Make the composite plot and tables for a given dataframe/set of phenotypes
# @param df A dataframe where rows are scans and columns are phenotypes and features
# @param phenotypes A vector of phenotypes to analyze
# @param colNames A vector of cleaned names to use for the table columns
# @param rowNames A vector of cleaned names to use for the table rows
# @param title A string to differentiate the figures for this dataframe vs. other dataframes
generateHemispherePlots <- function(df, phenotypes, phenotypesNoWS, title){
  # Initialize variables
  # Set up empty variables to generate table later
  ageAtPeak <- c()

  modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
                           sex='M', 
                           top_scan_reason_factors='headaches')
  
  # for phenotype in phenotypes...
  for (phenotype in phenotypesNoWS){
    # Generate GAMs using the formula that incorporates scan reason
    gamm <- createScanReasonGamm(df, phenotype)

    # Predict on the GAMs for the actual data
    gammPreds <- predict_gam(gamm$gam, values = modelFixedValues)

    # Get the age at peak
    ageAtPeak <- append(ageAtPeak, c(getAgeAtPeak(gammPreds)))
  }
  
  # Convert the lists into a dataframe
  results = as.data.frame(cbind(region=phenotypes, 
                                feat=phenotypesNoWS,
                                em=as.numeric(ageAtPeak)),
                       stringsAsFactors=F)
  
  print(results)
  
  # Plot the dataframe using ggseg
  p <- results %>% 
    ggseg(mapping=aes(fill=as.numeric(em)),
          hemisphere='left') +
    labs(title = "Age at Peak (years)", legend="Age") +
    theme(axis.title = element_blank()) +
    scale_fill_gradient(low = "blue", high = "red", na.value = NA) 
    
  grid.arrange(p)
  return(results)
}

regionalPeaks <- generateHemispherePlots(greyVolDf, regionalPhenotypes, parsedRegionalPhenotypes, '')


## Step 11: Compare CLIP age at peak to Lifespan age at peak -------------------
# Read the ages from lifespan
lifespanDf <- read.csv('/Users/youngjm/Data/lifespan_growth_charts/Lifespan_Data_Peaks_Table_2_2.csv')
# names(lifespanDf)[names(lifespanDf) == "feat"] <- "region"

regionalPeaks <- merge(regionalPeaks, lifespanDf, by='feat')
regionalPeaks$peakDiff <- as.numeric(regionalPeaks$Peak) - as.numeric(regionalPeaks$em)


regionalPeaks %>% 
  ggseg(mapping=aes(fill=Peak),
        hemisphere='left') +
  labs(title = "Lifespan Age at Peak (years)") +
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) 


# p <- lifespanDf %>%
#   ggseg(mapping=aes(fill=as.numeric(peakDiff)),
#         hemisphere='left') +
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
#   labs(title="Regional Peak Age\nDifference from Lifespan")
# 
# grid.arrange(p)
# 
ggplot(data=regionalPeaks, aes(x=region, y=peakDiff, fill=region)) +
  geom_bar(stat = 'identity', position = 'identity') +
  coord_flip() +
  theme_minimal() +
  guides(fill = 'none') +
  xlab('') + ylab('Difference in Age at Peak from Lifespan')

# Calculate the correlation between Lifespan and CLIP age at peak
tmp <- regionalPeaks[regionalPeaks$feat != 'medialorbitofrontal',]
ggplot(data=tmp, aes(color=as.factor(feat)))+
  geom_point(aes(x=as.numeric(em), y=Peak)) +
  geom_abline(slope = 1) +
  theme_minimal() +
  guides(fill = 'none') +
  # expand_limits(x=c(0, max(regionalPeaks$em)), y=c(0, max(regionalPeaks$Peak))) +
  xlab('Age at Peak CLIP (years)') +
  ylab('Age at Peak Lifespan (years)') +
  title('Age at Peak: Lifespan vs. CLIP')

cor(tmp$Peak, as.numeric(tmp$em))

## Step 12: Load SynthSeg metrics and generate GAMs ----------------------------

synthsegFn <- '/Users/youngjm/Data/clip/images/derivatives/metrics_synthsegplus.csv'
synthsegData <- read.csv(synthsegFn)

# Drop any rows where neurofibromatosis is in the scan_reason_categories column
synthsegData <- synthsegData[!grepl("neurofibromatosis", synthsegData$scan_reason_primary), ]
synthsegData <- synthsegData[!grepl("neurofibromatosis", synthsegData$scan_reason_categories), ]
synthsegData <- synthsegData[!grepl("true", synthsegData$confirm_neurofibromatosis), ]

# Make an age in years column
synthsegData$age_in_years <- synthsegData$age_in_days/365.25

# Add new column: going to sum TotalGrayVol + CerebralWhiteMatterVol + VentricleVolume + SubCortGrayVol
#analysisDf$TotalBrainVol <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol + analysisDf$VentricleVolume + analysisDf$SubCortGrayVol

# Some of the columns in this df should be factors
toFactor <- c('sex', 'Processing', 'MagneticFieldStrength', 'scanner_id', 
              'scan_reason_primary')
synthsegData[toFactor] <- lapply(synthsegData[toFactor], factor)
synthsegData <- synthsegData[(synthsegData$rawdata_image_grade >= 1),]
synthsegData <- addPrimaryScanReasonCol(synthsegData)


phenotypes <- c('SS_TotalBrainVolume', 'SS_TotalGrayVolume', 'SS_CerebralWhiteMatterVolume', 
                'SS_VentricleVolume', 'SS_SubcorticalGrayVolume')

# Build a table with all of the global phenotype p-vals parameters
tableColumns <- c("Total Brain\nVolume", "Total Gray\nVolume", "Cerebral White\nMatter Volume",
                  "Ventricle\nVolume", "Subcortical Gray\nVolume")
tableRows <- c("Surface Holes", "Scan Reason:\nDevelopmental Disorder",
               "Scan Reason:\nEye/Vision Finding", "Scan Reason:\nHeadaches",
               "Scan Reason:\nNon-Brain Lesion", "Scan Reason:\nSeizures",
               "Sex = Male", "Age")

rawAgeAtPeak <- generatePlotsAndTablesForDataset(synthsegData, phenotypes, tableColumns, tableRows, 'SynthSeg Phenotypes')

##### BONUS: Graphs for Aaron's grant ------------------------------------------
# Run steps 1-5, approximately


# Set up empty variables to generate table later
ageAtPeak <- list()
predictionPlots <- list()
tableValues <- list()
tableBetas <- list()

phenotypeTitles <- c("Total Brain Volume",
                     "Total Gray Volume",
                     "Cerebral White Matter Volume",
                     "Ventricle Volume",
                     "Subcortical Gray Volume",
                     "Cortical Surface Area",
                     "Mean Cortical Thickness")

modelFixedValues <- list(SurfaceHoles = mean(analysisDf$SurfaceHoles), 
                         sex='M', 
                         top_scan_reason_factors='headaches')

# for phenotype in phenotypes...
for (i in 1:length(phenotypes)){
  phenotype <- phenotypes[[i]]
  # Generate GAMs using the formula that incorporates scan reason
  gammScanReason <- createScanReasonGamm(analysisDf, phenotype)
  
  # Predict on the GAMs for the actual data
  gammScanReasonPreds <- predict_gam(gammScanReason$gam, values = modelFixedValues)
  
  # Get the age at peak
  ageAtPeak[[phenotype]] <- getAgeAtPeak(gammScanReasonPreds)
  
  # Generate scatter plots with confidence intervals
  plotScatterScanReason <- generatePlotScatterWithCI(gammScanReasonPreds, analysisDf, analysisDf[ , phenotype], phenotypeTitles[[i]])
  predictionPlots[[phenotype]] <- plotScatterScanReason
  
  # Add another plot with just the confidence intervals
  plotCI <- generateDiagnosisPlotCI(gammScanReasonPreds, analysisDf, analysisDf[ , phenotype], phenotypeTitles[[i]])
  predictionPlots[[paste(phenotype, 'ci', sep='_')]] <- plotCI
  
  tParam <- generateScanReasonsParametricTable(gammScanReason$gam)
  tLin <- generateScanReasonsLinearTable(gammScanReason$gam)
  tableValues[[phenotype]] <- append(tLin, list(Age = tParam))
  tableBetas[[phenotype]] <- getBetasFromGamSummary(gammScanReason$gam)
}

# Plot all scatter/CI plots in 1 figure
patchwork <- wrap_plots(predictionPlots, nrow=2, guides="collect", byrow = FALSE)
print(patchwork + plot_annotation(title=paste("Lifespan Trajectories of Global Phenotypes")))

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