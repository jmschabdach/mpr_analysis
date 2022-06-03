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
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
# globalPhenotypes <- c('TotalBrainVol', 'TotalGrayVol', 'CerebralWhiteMatterVol', 
#                       'VentricleVolume', 'SubCortGrayVol', 'CorticalSurfaceArea',
#                       'MeanCorticalThickness')

#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
#-------------------------------------------------------------------------------

## Demographic plot functions:
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

generateEulerQcDistributionPlot <- function(dataDf, sex, colorPalette){
  # Set the title string
  if (sex == 'M'){
    title = "Euler Number vs. QC Rating (Male)"
    dataDf <- dataDf[dataDf$sex == 'M', ]
  } else {
    title = "Euler Number vs. QC Rating (Female)"
    dataDf <- dataDf[dataDf$sex == 'F', ]
  }
  
  # Make plots for QC
  pQcSurfaceHoles <- dataDf %>%
    ggplot(aes(x=as.factor(rawdata_image_grade), y=SurfaceHoles)) +
    geom_violin(aes(fill=as.factor(rawdata_image_grade))) +
    geom_boxplot(width=0.1)+ 
    scale_fill_manual(values=qcColors) +
    labs(title = title,
         x = 'Image QC Rating Group',
         y = 'Number of Surface Holes (count)',
         fill="Image Quality") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  
  # Estimate the effect?
}

generateAgeQcDistributionPlot <- function(dataDf, sex, colorPalette){
  
  # Set the title string
  if (sex == 'M'){
    title = "Age vs. QC Rating (Male)"
    dataDf <- dataDf[dataDf$sex == 'M', ]
  } else {
    title = "Age vs. QC Rating (Female)"
    dataDf <- dataDf[dataDf$sex == 'F', ]
  }
  
  pQcAge <- dataDf %>%
    ggplot(aes(x=as.factor(rawdata_image_grade), y=age_in_years)) +
    geom_violin(aes(fill=as.factor(rawdata_image_grade))) +
    geom_boxplot(width=0.1)+ 
    scale_fill_manual(values=qcColors)+
    labs(title = title,
         x = 'Image QC Rating Group',
         y = 'Age at Scan (Days)',
         fill="Image Quality") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
}

##
# Given a dataframe and a column header (imaging phenotype), create a GAMM
# (GAMM = Generalized Additive Mixed Model)
# @param df A dataframe where the rows are scans and columns contain phenotypes
# @param measure Name of a column in the dataframe
# @return mixedModel The mixed model for the given formula and specificed measure
createMainGamm <- function(df, measure) {
  formula <- as.formula(paste(measure, "s(log(age_in_years), k=3) +
                      SurfaceHoles +
                      sex", sep="~"))
  
  # Make a basic linear model accounting for age and surface holes for the given data frame
  mixedModel <- gam(formula,
                    data = df,
                    gamma=1)
  # Return the model
  return(mixedModel)
}

createGammLinearAge <- function(df, measure) {
  formula <- as.formula(paste(measure, "log(age_in_years) +
                      SurfaceHoles +
                      # top_scan_reason_factors +
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

createGammScanReasonsAge <- function(df, measure) {
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

createGammNoSex <- function(df, measure) {
  formula <- as.formula(paste(measure, "s(log(age_in_years), fx=T) +
                      SurfaceHoles", sep="~"))
  
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
    geom_point(data=origData, alpha=0.5, aes(x=age_in_years, y=measure, color=sex)) +
    geom_smooth_ci() +
    geom_vline(xintercept=peakAge, linetype='dashed', color='#c23400') +
    scale_color_manual(values = cbbPalette, name = "Sex") +
    # scale_x_continuous(trans='log10') +
    theme(axis.title = element_blank(), plot.title=element_text(hjust=0.5)) +
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
  print(subTab)

  newRowNames <- c('Surface Holes',
                   # 'Scan Reason: Developmental Disorder',
                   # 'Scan Reason: Eye/Vision Finding',
                   # 'Scan Reason: Headaches',
                   # 'Scan Reason: Non-Brain Lesion',
                   # 'Scan Reason: Seizures',
                   'Sex')
  rows <- rownames(subTab)[2:length(rownames(subTab))]
  subTab <- subset(subTab, rownames(subTab) %in% rows)
  rownames(subTab) <- newRowNames

  # Rename and trim the columns
  colnames(subTab) <- c('beta', 'std err', 't-value', 'p-value')
  cols <- c('p-value')
  subTab <- subTab[, colnames(subTab) %in% cols]

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
                   # 'Scan Reason: Developmental Disorder',
                   # 'Scan Reason: Eye/Vision Finding',
                   # 'Scan Reason: Headaches',
                   # 'Scan Reason: Non-Brain Lesion',
                   # 'Scan Reason: Seizures',
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
  for (i in 1:length(phenotypes)){
    print(i)
    phenotype <- phenotypes[[i]]
    col <- colNames[i]
    print(col)
    # Generate GAMs using the formula that incorporates scan reason
    gammScanReason <- createMainGamm(df, phenotype)
    
    
    # Predict on the GAMs for the actual data
    # gammScanReasonPreds <- predict_gam(gammScanReason$gam, values = modelFixedValues)
    gammScanReasonPreds <- predict_gam(gammScanReason, values = modelFixedValues)
    
    
    # Get the age at peak
    ageAtPeak[[phenotype]] <- getAgeAtPeak(gammScanReasonPreds)
    
    # Generate scatter plots with confidence intervals
    plotScatterScanReason <- generatePlotScatterWithCI(gammScanReasonPreds, df, df[ , phenotype], col)
    predictionPlots[[phenotype]] <- plotScatterScanReason
    
    # tParam <- generateScanReasonsParametricTable(gammScanReason$gam)
    # tLin <- generateScanReasonsLinearTable(gammScanReason$gam)
    tParam <- generateScanReasonsParametricTable(gammScanReason)
    tLin <- generateScanReasonsLinearTable(gammScanReason)
    tableValues[[phenotype]] <- append(tLin, list(Age = tParam))
    
    # Get normalized betas
    normDf <- data.frame(df)
    normDf[[phenotype]] <- scale(normDf[[phenotype]])
    normDf$SurfaceHoles <- scale(normDf$SurfaceHoles)
    normDf$age_in_years <- scale(normDf$age_in_years)
    gammNorm <- createMainGamm(normDf, phenotype)
    # tableBetas[[phenotype]] <- getBetasFromGamSummary(gammNorm$gam)
    tableBetas[[phenotype]] <- getBetasFromGamSummary(gammNorm)
  }
  
  # Plot all scatter/CI plots in 1 figure
  patchwork <- wrap_plots(predictionPlots, nrow=1, guides="collect")
  print(patchwork + plot_annotation(title=paste("Trajectories of", title)))
  # Make table: p-values of factor/phenotype
  firstStep <- lapply(tableValues, unlist)
  pvals <- as.data.frame(firstStep, stringsAsFactors = F)
  colnames(pvals) <- colNames
  rownames(pvals) <- rowNames
  # t1 <- gridExtra::tableGrob(secondStep)
  # gridExtra::grid.arrange(top=paste("p-values of", title), t1)

  # Make table: betas
  firstStep <- lapply(tableBetas, unlist)
  betaVals <- as.data.frame(firstStep, stringsAsFactors = F)
  colnames(betaVals) <- tableColumns
  rownames(betaVals) <- tableRows[1:length(tableRows)-1]
  # t2 <- gridExtra::tableGrob(secondStep)
  # gridExtra::grid.arrange(top=paste("Betas of", title), t1)
  
  output <- list()
  output$predictionPlots <- predictionPlots
  output$pvalsTable <- pvals
  output$betasTable <- betaVals
  output$ageAtPeak <- ageAtPeak
  
  return(output)
}

##
# Make the composite plot and tables for a given dataframe/set of phenotypes
# @param df A dataframe where rows are scans and columns are phenotypes and features
# @param phenotypes A vector of phenotypes to analyze
# @param colNames A vector of cleaned names to use for the table columns
# @param rowNames A vector of cleaned names to use for the table rows
# @param title A string to differentiate the figures for this dataframe vs. other dataframes
generateCentilePlots <- function(df, phenotypes, cols, title){
  # Initialize variables
  centilePlots <- list()
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # for phenotype in phenotypes...
  for (i in 1:length(phenotypes)){
    phenotype <- phenotypes[[i]]
    # Generate scatter plots with confidence intervals
    plot01 <- df %>%
      ggplot(aes(x=age_in_years, y=.data[[phenotype]]), log='x') +
      geom_point(data=df, alpha=0.5, aes(color=sex)) +
      scale_color_manual(values = cbbPalette, name = "Sex") +
      theme(axis.title = element_blank(), plot.title=element_text(hjust=0.5)) +
      labs(title = cols[[i]])
    
    centilePlots[[phenotype]] <- plot01
  }
  
  # Plot all scatter/CI plots in 1 figure
  patchwork <- wrap_plots(centilePlots, nrow=1, guides="collect")
  print(patchwork + plot_annotation(title=paste("Centiles of", title)))
  
  return(centilePlots)
}


generatePlotsAndTablesLinearAge <- function(df, phenotypes, colNames, rowNames, title){
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
    gammScanReason <- createGammLinearAge(df, phenotype)
    
    # Predict on the GAMs for the actual data
    gammScanReasonPreds <- predict_gam(gammScanReason$gam, values = modelFixedValues)
    
    # Get the age at peak
    ageAtPeak[[phenotype]] <- getAgeAtPeak(gammScanReasonPreds)
    
    # Generate scatter plots with confidence intervals
    plotScatterScanReason <- generatePlotScatterWithCI(gammScanReasonPreds, df, df[ , phenotype], phenotype)
    predictionPlots[[phenotype]] <- plotScatterScanReason
    
    tLin <- signif(summary(gammScanReason$gam)$p.table, 3)
    colnames(tLin) <- c('beta', 'std err', 't-value', 'p-value')
    tableValues[[phenotype]] <- tLin
  }
  
  # Plot all scatter/CI plots in 1 figure
  patchwork <- wrap_plots(predictionPlots, nrow=1, guides="collect")
  print(patchwork + plot_annotation(title=paste("Trajectories of", title)))
  
  return(tableValues)
}


#-------------------------------------------------------------------------------
# Loading and Prepping Data
#-------------------------------------------------------------------------------

## Step 1: load data and prep it ------------------------------------------------

# Load the master file containing subject demographics, imaging phenotypes, etc.
rawFn <- '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-03_analysis_features.csv'
fs6Fn <- '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_structural_stats.csv'

# inFn <- '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_7.1.1_structural_stats.csv'
masterDf <- read.csv(fs6Fn)

# Drop any rows where neurofibromatosis is in the scan_reason_categories column
analysisDf <- masterDf[!grepl("neurofibromatosis", masterDf$scan_reason_primary), ]
analysisDf <- analysisDf[!grepl("neurofibromatosis", analysisDf$scan_reason_categories), ]
# Need to convert missing values in "confirm_neurofibromatosis" to FALSE
analysisDf <- analysisDf %>%
  mutate(confirm_neurofibromatosis = case_when(
    grepl('1', confirm_neurofibromatosis, fixed=TRUE) ~ TRUE,
    grepl('0', confirm_neurofibromatosis, fixed=TRUE) ~ FALSE,
    is.na(confirm_neurofibromatosis) ~ FALSE
  ))
analysisDf <- analysisDf[(analysisDf$confirm_neurofibromatosis == FALSE), ]

# Make an age in years column
analysisDf$age_in_years <- analysisDf$age_at_scan_days/365.25

# Add new column: going to sum TotalGrayVol + CerebralWhiteMatterVol + VentricleVolume + SubCortGrayVol
analysisDf$TotalBrainVol <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol + analysisDf$VentricleVolume + analysisDf$SubCortGrayVol

# Add a column for processing
analysisDf$Processing <- 'FS6'

# Some of the columns in this df should be factors
toFactor <- c('sex', 'Processing', 'MagneticFieldStrength', 'scanner_id', 
              'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

# Drop any 1.5T scans
analysisDf <- analysisDf[analysisDf$MagneticFieldStrength != "1.5",]

# Drop any scans with ratings less than 0
analysisDf <- analysisDf[analysisDf$rawdata_image_grade >= 0, ]

# Only one scan per subject
# Sort the dataframe by patient_id and scanner_id
analysisDf <- analysisDf[ with(analysisDf, order(analysisDf$patient_id, analysisDf$scan_id)), ]
analysisDf <- analysisDf[!duplicated(analysisDf$patient_id), ]
analysisDf$patient_id <- droplevels(as.factor(analysisDf$patient_id))

# Add a column for TCV (Total Cerebrum Volume)
analysisDf$TCV <- analysisDf$TotalGrayVol + analysisDf$CerebralWhiteMatterVol

# Drop any scans with NAs
analysisDf <- analysisDf[complete.cases(analysisDf), ]
write.csv(analysisDf, '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-05-26_highres_nocontrast_singlescanpersubject.csv')

## Step 2: Generate basic demographic plots ------------------------------------
# These plots should be consistent across all analyses

qcColors <- c("#999999", "#0072B2", "#56B4E9")

# Make a plot for the distribution of ages with coloring for sex 
pAgesMale <- generateAgeDistributionPlot(analysisDf, 'M', qcColors)
pAgesFemale <- generateAgeDistributionPlot(analysisDf, 'F', qcColors)

pScannersMale <- generateScannerDistributionPlot(analysisDf, 'M', qcColors)
pScannersFemale <- generateScannerDistributionPlot(analysisDf, 'F', qcColors)

# Make plots for QC
pEulerQcMale <- generateEulerQcDistributionPlot(analysisDf, 'M', qcColors)
pEulerQcFemale <- generateEulerQcDistributionPlot(analysisDf, 'F', qcColors)

pAgeQcMale <- generateAgeQcDistributionPlot(analysisDf, 'M', qcColors)
pAgeQcFemale <- generateAgeQcDistributionPlot(analysisDf, 'F', qcColors)

ttest2v1M <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 1, ]$age_in_years)
ttest1v0M <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 1, ]$age_in_years, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v0M <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v1F <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 1, ]$age_in_years)
ttest1v0F <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 1, ]$age_in_years, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v0F <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v1All <- t.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$rawdata_image_grade == 1, ]$age_in_years)
ttest1v0All <- t.test(analysisDf[analysisDf$rawdata_image_grade == 1, ]$age_in_years, analysisDf[analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v0All <- t.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$rawdata_image_grade == 0, ]$age_in_years)


qcAgeTableM <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                          pvalues = c(ttest1v0M$p.value, ttest2v0M$p.value, ttest2v1M$p.value),
                          ciLower = c(ttest1v0M$conf.int[1], ttest2v0M$conf.int[1], ttest2v1M$conf.int[1]),
                          ciUpper = c(ttest1v0M$conf.int[2], ttest2v0M$conf.int[2], ttest2v1M$conf.int[2]))

qcAgeTableF <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                          pvalues = c(ttest1v0F$p.value, ttest2v0F$p.value, ttest2v1F$p.value),
                          ciLower = c(ttest1v0F$conf.int[1], ttest2v0F$conf.int[1], ttest2v1F$conf.int[1]),
                          ciUpper = c(ttest1v0F$conf.int[2], ttest2v0F$conf.int[2], ttest2v1F$conf.int[2]))

qcAgeTableAll <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                          pvalues = c(ttest1v0All$p.value, ttest2v0All$p.value, ttest2v1All$p.value),
                          ciLower = c(ttest1v0All$conf.int[1], ttest2v0All$conf.int[1], ttest2v1All$conf.int[1]),
                          ciUpper = c(ttest1v0All$conf.int[2], ttest2v0All$conf.int[2], ttest2v1All$conf.int[2]))

# Test for Euler and QC
ttest2v1EulerM <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles)
ttest1v0EulerM <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v0EulerM <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v1EulerF <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles)
ttest1v0EulerF <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v0EulerF <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v1EulerAll <- t.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles)
ttest1v0EulerAll <- t.test(analysisDf[analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles, analysisDf[analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v0EulerAll <- t.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)

qcAgeTableEulerM <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                          pvalues = c(ttest1v0EulerM$p.value, ttest2v0EulerM$p.value, ttest2v1EulerM$p.value),
                          ciLower = c(ttest1v0EulerM$conf.int[1], ttest2v0EulerM$conf.int[1], ttest2v1EulerM$conf.int[1]),
                          ciUpper = c(ttest1v0EulerM$conf.int[2], ttest2v0EulerM$conf.int[2], ttest2v1EulerM$conf.int[2]))

qcAgeTableEulerF <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                          pvalues = c(ttest1v0EulerF$p.value, ttest2v0EulerF$p.value, ttest2v1EulerF$p.value),
                          ciLower = c(ttest1v0EulerF$conf.int[1], ttest2v0EulerF$conf.int[1], ttest2v1EulerF$conf.int[1]),
                          ciUpper = c(ttest1v0EulerF$conf.int[2], ttest2v0EulerF$conf.int[2], ttest2v1EulerF$conf.int[2]))

qcAgeTableEulerAll <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                            pvalues = c(ttest1v0EulerAll$p.value, ttest2v0EulerAll$p.value, ttest2v1EulerAll$p.value),
                            ciLower = c(ttest1v0EulerAll$conf.int[1], ttest2v0EulerAll$conf.int[1], ttest2v1EulerAll$conf.int[1]),
                            ciUpper = c(ttest1v0EulerAll$conf.int[2], ttest2v0EulerAll$conf.int[2], ttest2v1EulerAll$conf.int[2]))



# Make a table containing the number of scans with each rating divided by sex
qcCountTable <- data.frame(rating = c("2", "1", "0"),
                           Male = c(dim(analysisDf[analysisDf$sex == "M" & analysisDf$rawdata_image_grade == 2, ])[1],
                                    dim(analysisDf[analysisDf$sex == "M" & analysisDf$rawdata_image_grade == 1, ])[1],
                                    dim(analysisDf[analysisDf$sex == "M" & analysisDf$rawdata_image_grade == 0, ])[1]),
                           Female = c(dim(analysisDf[analysisDf$sex == "F" & analysisDf$rawdata_image_grade == 2, ])[1],
                                      dim(analysisDf[analysisDf$sex == "F" & analysisDf$rawdata_image_grade == 1, ])[1],
                                      dim(analysisDf[analysisDf$sex == "F" & analysisDf$rawdata_image_grade == 0, ])[1]))



grob <- patchworkGrob(pAgesMale + pAgesFemale + 
                        pScannersMale + pScannersFemale +
                        pAgeQcMale + pAgeQcFemale +
                        pEulerQcMale + pEulerQcFemale +
                        plot_layout(guides="collect", ncol = 2))
gridExtra::grid.arrange(grob)

# ^^^ SAVE THIS FIGURE


# Plot the distribution of reasons for a scan
analysisDf <- addPrimaryScanReasonCol(analysisDf)
reasonsTable <- summary(analysisDf$top_scan_reason_factors)
pieLabels <- paste(names(reasonsTable), reasonsTable, sep='\n')
pie(reasonsTable, pieLabels, col=cbbPalette, mai=c(0,0,0,0),
    main="Top Scan Reasons for Clinical Controls\n(with sample sizes)")

## Drop the rows where rawdata_image_grade <= 1 --------------------------------
highQDf <- analysisDf[(analysisDf$rawdata_image_grade >= 1),]
superHighQDf <- analysisDf[analysisDf$rawdata_image_grade >1 , ]

## Step 3: ComBat the data/Load ComBat data ------------------------------------
prepForCombat <- function(df, fnBase){
  # Move metadata to the front of the dataframe
  cols <- colnames(df)
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",                
                "MagneticFieldStrength", "scanner_id", "confirm_neurofibromatosis", 
                "rawdata_image_grade", 'Processing', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories", 'SurfaceHoles')
  phenoCols <- setdiff(cols, metaCols)
  globalPhenoCols <- c('TotalGrayVol', 'CerebralWhiteMatterVol', 'SubCortGrayVol',
                       'eTIV', 'VentricleVolume', 'CorticalSurfaceArea', 
                       'MeanCorticalThickness', 'TCV')
  nonGlobalCols <- setdiff(phenoCols, globalPhenoCols)
  
  regionalPhenoCols <- c()
  for (c in nonGlobalCols){
    if ((grepl('lh_', c) | grepl('rh_', c)) & grepl('_grayVol', c)){
      regionalPhenoCols <- append(regionalPhenoCols, c)
    }
  }
  
  newCols <- c(metaCols, globalPhenoCols, regionalPhenoCols)
  df <- df[, newCols]
  
  # Drop any scans with NAs
  df <- df[complete.cases(df), ]
  
  # Identify the scan_id + phenotypes to combat
  toCombat <- df[, c('scan_id', globalPhenoCols, regionalPhenoCols)]
  
  # Identify covariates to harmonize on/protect
  covars <- data.frame(SITE = df$scanner_id, # <-- harmonize on the first column
                       log_age_in_years = log(df$age_in_years), # <-- protect all following columns
                       sex = df$sex,
                       reason = df$top_scan_reason_factors)
  
  # Save the data and covars dataframes to csvs
  write.csv(toCombat,paste0(fnBase, "_phenotypes_toCombat.csv"), row.names = FALSE)
  write.csv(covars,paste0(fnBase, "_covariates_toCombat.csv"), row.names = FALSE)
}

prepForCombat(analysisDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc0-2")
prepForCombat(highQDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc1-2")
prepForCombat(superHighQDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc2")

# STOP HERE AND RUN COMBAT VIA PYTHON

loadCombattedData <- function(df, fn){
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",                
                "MagneticFieldStrength", "scanner_id", "confirm_neurofibromatosis", 
                "rawdata_image_grade", 'Processing', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories", 'SurfaceHoles')
  combattedDf <- read.csv(fn)
  # Drop scan_ids column (duplicate)
  combattedDf <- combattedDf[ , -which(names(combattedDf) %in% c("scan_id"))]
  # Add the metadata back in
  metadataDf <- df[, metaCols]
  combattedDf <- cbind(metadataDf, combattedDf)
  return(combattedDf)
}

combattedDf <- loadCombattedData(analysisDf, '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc0-2_combatted.csv')
combattedHighQDf <- loadCombattedData(highQDf, '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc1-2_combatted.csv')
combattedSuperHighQDf <- loadCombattedData(superHighQDf, '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc2_combatted.csv')

combattedDf <- combattedDf[complete.cases(combattedDf), ]
combattedHighQDf <- combattedHighQDf[complete.cases(combattedHighQDf), ]
combattedSuperHighQDf <- combattedSuperHighQDf[complete.cases(combattedSuperHighQDf), ]

## Step 4: Prep data for centilization ---------

prepForCentilizing <- function(df, fn){
  toCentilizeCols <- c('patient_id', 'age_in_years', 'sex', 'TotalGrayVol', 
          'CerebralWhiteMatterVol', 'SubCortGrayVol', 'VentricleVolume',
          'CorticalSurfaceArea', 'MeanCorticalThickness', 'TCV',
          'SurfaceHoles', 'scanner_id', 'top_scan_reason_factors')
  newNames <- c('participant', 'Age', 'sex', 'GMV', 'WMV', 'sGMV', 'Ventricles', 
                'SA', 'CT', 'TCV',
                'SurfaceHoles', 'scanner_id', 'top_scan_reason_factors')
  toCentilizeDf <- df[, toCentilizeCols]
  toCentilizeDf <- setnames(toCentilizeDf, toCentilizeCols, newNames)
  
  toCentilizeDf <- toCentilizeDf %>%
    mutate(sex = case_when(
      grepl('F', sex, fixed=TRUE) ~ "Female",
      grepl('M', sex, fixed=TRUE) ~ "Male"
    ))
  
  toCentilizeDf$age_days <- df$age_at_scan_days + 280 # add gestational age
  toCentilizeDf$study <-"JMY_CONT_FS6"
  toCentilizeDf$fs_version <- "FS6_T1"
  toCentilizeDf$country <- "USA"
  toCentilizeDf$run <- 1
  toCentilizeDf$session <- 1
  toCentilizeDf$dx <- "CN"
  reorderCols <- c('participant', 'Age', 'age_days', 'sex', 'study', 'fs_version',
            'scanner_id', 'country', 'run', 'session', 'dx', 'top_scan_reason_factors',
            'SurfaceHoles',
            'GMV', 'WMV', 'sGMV', 'Ventricles', 'SA', 'CT', 'TCV')
  toCentilizeDf <- toCentilizeDf[, reorderCols]
  toCentilizeDf$INDEX.TYPE <- "NA"
  toCentilizeDf$INDEX.OB <- "NA"
  
  write.csv(toCentilizeDf, fn)
  return(toCentilizeDf)
}

preCenDf <- prepForCentilizing(analysisDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc0-2_raw_to_centilize.csv")
preCenHighQDf <- prepForCentilizing(highQDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc1-2_raw_to_centilize.csv")
preCenSuperHighQDf <- prepForCentilizing(superHighQDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc2_raw_to_centilize.csv")

comPreCenDf <- prepForCentilizing(combattedDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc0-2_combatted_to_centilize.csv")
comPreCenHighQDf <- prepForCentilizing(combattedHighQDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc1-2_combatted_to_centilize.csv")
comPreCenSuperHighQDf <- prepForCentilizing(combattedSuperHighQDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc2_combatted_to_centilize.csv")


## Step 5: Load centilized phenotypes ---------------------------------

loadCentilizedFeature <- function(centileFn, phenotype){
  centileDf <- read.csv(centileFn)
  print(head(centileDf))

  c1 <- paste0(phenotype, "Transformed.normalised")
  c2 <- paste0(phenotype, "Transformed.q.wre")
  print(c(c1, c2))
  idx1 <- grep(c1, colnames(centileDf))[1]
  idx2 <- grep(c2, colnames(centileDf))[1]
  col <- c(colnames(centileDf)[idx1], colnames(centileDf[idx2]))

  return(centileDf[col])
}

loadCentilizedData <- function(centileDf, fnsBase){
  centileCols <- c('GMV', 'WMV', 'sGMV', 
                   'Ventricles', 'totalSA2', 'meanCT2',
                   'TCV')
  
  centileFns <- c(paste0(fnsBase, '_GMV.csv'),
                  paste0(fnsBase, '_WMV.csv'),
                  paste0(fnsBase, '_sGMV.csv'),
                  paste0(fnsBase, '_Ventricles.csv'),
                  paste0(fnsBase, '_SA.csv'),
                  paste0(fnsBase, '_CT.csv'),
                  paste0(fnsBase, '_TCV.csv'))
  
  for (i in (1:length(centileFns))){
    # centileDf <- cbind(centileDf, loadCentilizedFeature(centileFns[[i]], paste0(centileCols[[i]], 'Transformed.q.wre')))
    tmpDf <- loadCentilizedFeature(centileFns[[i]], centileCols[[i]])
    # print(tmpDf)
    centileDf <- cbind(centileDf, tmpDf)
    
    names(centileDf)[names(centileDf) == 'centile'] <- paste0(centileCols[[i]], '_centile')
  }
  
  setnames(centileDf, old=c('participant', 'Age'), new=c('patient_id', 'age_in_years'))
  centileDf <- centileDf %>%
    mutate(sex = case_when(
      grepl('Female', sex, fixed=TRUE) ~ "F",
      grepl('Male', sex, fixed=TRUE) ~ "M"
    ))
  
  # Reorder scan reason factor levels
  centileDf$top_scan_reason_factors <- as.factor(centileDf$top_scan_reason_factors)
  centileDf$top_scan_reason_factors <- relevel(centileDf$top_scan_reason_factors, "other")
  
  return(centileDf)
}

ccDf <- loadCentilizedData(comPreCenDf, '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc0-2_combatted_centilized')
# ccHighQDf <- loadCentilizedData(comPreCenHighQDf, '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc1-2_combatted_centilized')
# ccSuperHighQDf <- loadCentilizedData(comPreCenSuperHighQDf, '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/qc2_combatted_centilized')

## Step 6: Generate GAMs for each phenotype of interest ------------------------

# Build a table with all of the global phenotype p-vals parameters
tableColumns <- c("Total Gray\nVolume", "Cerebral White\nMatter Volume", 
                  "Subcortical Gray\nVolume", "Ventricle\nVolume", 
                  "Total Cortical\nSurface Area", "Mean Cortical\nThickness", 
                  "Total Cerebrum\nVolume") 

tableRows <- c("Surface Holes", 
               "Sex = Male", 
               "Age")

miniGlobalCols <- c("TotalGrayVol", "CerebralWhiteMatterVol", "SubCortGrayVol",
                    "VentricleVolume", "CorticalSurfaceArea", 
                    "MeanCorticalThickness", "TCV")

rawPlotsTables <- generatePlotsAndTablesForDataset(highQDf, miniGlobalCols, tableColumns, tableRows, 'Raw Phenotypes')
combatPlotsTables <- generatePlotsAndTablesForDataset(combattedDf, miniGlobalCols, tableColumns, tableRows, 'Combatted Phenotypes')
combatHighQPlotsTables <- generatePlotsAndTablesForDataset(combattedHighQDf, miniGlobalCols, tableColumns, tableRows, 'Combatted Phenotypes')
combatSuperHighQPlotsTables <- generatePlotsAndTablesForDataset(combattedSuperHighQDf, miniGlobalCols, tableColumns, tableRows, 'Combatted Phenotypes')

centileCols <- c('GMVTransformed.normalised', 'WMVTransformed.normalised', 
                 'sGMVTransformed.normalised', 'VentriclesTransformed.normalised', 
                 "totalSA2Transformed.normalised", 
                  "meanCT2Transformed.normalised", "TCVTransformed.normalised")
# centiledRawPlotsTables <- generatePlotsAndTablesForDataset(centiledRawDf, centileCols, tableColumns, tableRows, 'Centilized Raw Phenotypes')
centiledCombattedPlotsTables <- generatePlotsAndTablesForDataset(ccDf, centileCols, tableColumns, tableRows, 'Centilized Combatted Phenotypes')
# centiledCombattedHighQPlotsTables <- generatePlotsAndTablesForDataset(ccHighQDf, centileCols, tableColumns, tableRows, 'Centilized Combatted Phenotypes')
# centiledCombattedSuperHighQPlotsTables <- generatePlotsAndTablesForDataset(ccSuperHighQDf, centileCols, tableColumns, tableRows, 'Centilized Combatted Phenotypes')

## Results 1: Plot the centiles ---------------------------------------------
centileList <- c("GMVTransformed.q.wre", "WMVTransformed.q.wre",
                 "sGMVTransformed.q.wre", "VentriclesTransformed.q.wre",
                 "totalSA2Transformed.q.wre", "meanCT2Transformed.q.wre",
                 "TCVTransformed.q.wre")
centileColTitles <- c("Total Gray\nVolume Centiles", 
                  "Cerebral White Matter\nVolume Centiles", 
                  "Subcortical Gray\nVolume Centiles", 
                  "Ventricle Volume\nCentiles", 
                  "Total Cortical Surface\nArea Centiles", 
                  "Mean Cortical\nThickness Centiles", 
                  "Total Cerebrum\nVolume Centiles")
ccCentilePlots <- generateCentilePlots(ccDf, centileList, centileColTitles, "Post Combat Centiles for QC 0-2")
# ccCentileHighQPlots <- generateCentilePlots(ccHighQDf, centileList, centileColTitles, "Post Combat Centiles for QC 1-2")
# ccCentileSuperHighQPlots <- generateCentilePlots(ccSuperHighQDf, centileList, centileColTitles, "Post Combat Centiles for QC 2")

# combatCentileHighQPlotsTables <- generatePlotsAndTablesForDataset(ccHighQDf, centileList, centileColTitles, tableRows, 'Centiles of Combatted Phenotypes')


## Results 2: Make tables ---------------------------------------------------------

# # Age at peak table
# rawAgeAtPeaks <- as.data.frame(rawPlotsTables$ageAtPeak, row.names = c("Raw Phenotypes"))
# combattedAgeAtPeaks <- as.data.frame(combatPlotsTables$ageAtPeak, row.names = c("Combatted Phenotypes"))
# centilizedRawAgeAtPeaks <- as.data.frame(centiledRawPlotsTables$ageAtPeak, row.names = c("Adjusted Raw Phenotypes"))
# centilizedCombattedAgeAtPeaks <- as.data.frame(centiledCombattedPlotsTables$ageAtPeak, row.names = c("Adjusted Combatted Phenotypes"))
# colnames(centilizedRawAgeAtPeaks) <- colnames(rawAgeAtPeaks)
# colnames(centilizedCombattedAgeAtPeaks) <- colnames(rawAgeAtPeaks)
# 
# agesAtPeakDf <- rbind(rawAgeAtPeaks, combattedAgeAtPeaks,
#                       centilizedRawAgeAtPeaks, centilizedCombattedAgeAtPeaks)

# Make stats table with p values
rawPvals <- as.data.frame(rawPlotsTables$pvalsTable) #, row.names = c("Raw Phenotypes"))
combattedPvals <- as.data.frame(combatHighQPlotsTables$pvalsTable) #, row.names = c("Combatted Phenotypes"))
centilePvals <- as.data.frame(combatCentileHighQPlotsTables$pvalsTable)
# centilizedPvals <- as.data.frame(centiledRawPlotsTables$pvalsTable) #, row.names = c("Adjusted Raw Phenotypes"))
# centilizedCombattedPvals <- as.data.frame(centiledCombattedPlotsTables$pvalsTable) #, row.names = c("Adjusted Combatted Phenotypes"))

rawPvals$Category <- "Raw Phenotypes"
combattedPvals$Category <- "Combatted Phenotypes"
centilePvals$Category <- "Centiles"
# centilizedPvals$Category <- "Centilized Phenotypes"
# centilizedCombattedPvals$Category <- "Centilized Combatted Phenotypes"

pvalsDf <- rbind(rawPvals, combattedPvals) #, centilizedPvals, centilizedCombattedPvals)
options(scipen = 0)
pvalsDf
centilePvals
# write.csv(pvalsDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/2022-05-16_global_pvals.csv")

# Make stats table with p values
rawBetas <- as.data.frame(rawPlotsTables$betasTable) #, row.names = c("Raw Phenotypes"))
combattedBetas <- as.data.frame(combatHighQPlotsTables$betasTable) #, row.names = c("Combatted Phenotypes"))
centileBetas <- as.data.frame(combatCentileHighQPlotsTables$betasTable)
# centilizedBetas <- as.data.frame(centiledRawPlotsTables$betasTable) #, row.names = c("Adjusted Raw Phenotypes"))
# centilizedCombattedBetas <- as.data.frame(centiledCombattedPlotsTables$betasTable) #, row.names = c("Adjusted Combatted Phenotypes"))

rawBetas$Category <- "Raw Phenotypes"
combattedBetas$Category <- "Combatted Phenotypes"
# centilizedBetas$Category <- "Centilized Phenotypes"
# centilizedCombattedBetas$Category <- "Centilized Combatted Phenotypes"

# betasDf <- rbind(rawBetas, combattedBetas, centilizedBetas, centilizedCombattedBetas)
options(scipen = 0)
rawBetas
combattedBetas
centileBetas
# betasDf
# write.csv(betasDf, "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/2022-05-16_global_betas.csv")

## Results 3: Combine global phenotype plots ------------------------------------------------------

# patchwork <- wrap_plots(c(#rawPlotsTables$predictionPlots,
#                           combatPlotsTables$predictionPlots,
#                           #centiledRawPlotsTables$predictionPlots,
#                           centiledCombattedPlotsTables$predictionPlots),
#                         nrow=2, guides="collect")
# print(patchwork)

patchwork <- wrap_plots(c(combatPlotsTables$predictionPlots,
                          combatHighQPlotsTables$predictionPlots,
                          combatSuperHighQPlotsTables$predictionPlots),
                        ncol=3, byrow = FALSE, guides="collect")
print(patchwork)

patchwork <- wrap_plots(c(combatHighQPlotsTables$predictionPlots,
                          ccCentileHighQPlots,
                          centiledCombattedHighQPlotsTables$predictionPlots),
                        ncol=3, byrow = FALSE, guides="collect")
print(patchwork)

## Step 7: Prepare regional phenotype data -----------------------

brainDf <- combattedDf
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


grayVolCols <- brainDfCols[endsWith(brainDfCols, '_grayVol')]
print(grayVolCols)
phenoCols <- c()
for (pheno in parsedRegionalPhenotypes) {
  cols <- grayVolCols[grepl(paste("_",pheno, sep=''), grayVolCols, fixed=TRUE)]
  phenoCols <- c(phenoCols, cols)
}
grayVolLhCols <- sort(phenoCols[startsWith(phenoCols, 'lh_')])
grayVolRhCols <- sort(phenoCols[startsWith(phenoCols, 'rh_')])

# Let's calculate the bilateral average of each phenotype
for (i in 1:length(parsedRegionalPhenotypes)){
  brainDf[[parsedRegionalPhenotypes[[i]]]] <- (brainDf[[grayVolLhCols[[i]]]] + brainDf[[grayVolRhCols[[i]]]])/2
}

globalPhenoCols <- c('TotalGrayVol', 'CerebralWhiteMatterVol', 'SubCortGrayVol',
                     'eTIV', 'VentricleVolume', 'CorticalSurfaceArea', 
                     'MeanCorticalThickness', 'TCV')
regionalMetaCols <- setdiff(brainDfCols, c(globalPhenoCols, grayVolCols))
greyVolDf <- brainDf[append(regionalMetaCols, parsedRegionalPhenotypes)]

## Results 4: Make regional phenotype plots for age at peak --------------------
##
# Make the composite plot and tables for a given dataframe/set of phenotypes
# @param df A dataframe where rows are scans and columns are phenotypes and features
# @param phenotypes A vector of phenotypes to analyze
# @param colNames A vector of cleaned names to use for the table columns
# @param rowNames A vector of cleaned names to use for the table rows
# @param title A string to differentiate the figures for this dataframe vs. other dataframes
generateHemispherePlots <- function(df, phenotypes, phenotypesNoWS, title, sex=""){
  # Initialize variables
  # Set up empty variables to generate table later
  ageAtPeak <- c()
  
  # for phenotype in phenotypes...
  for (phenotype in phenotypesNoWS){
    if (sex == "M") {
      print("M")
      modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
                               sex="M",
                               top_scan_reason_factors='headaches')
      gamm <- createGammNoSex(df[df$sex == "M", ], phenotype)
      gammPreds <- predict_gam(gamm$gam, values = modelFixedValues)
    } else if (sex == "F"){
      print("F")
      modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
                               sex="F",
                               top_scan_reason_factors='headaches')
      gamm <- createGammNoSex(df[df$sex == "F", ], phenotype)
      gammPreds <- predict_gam(gamm$gam, values = modelFixedValues)
    } else {
      modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
                               sex='M',
                               top_scan_reason_factors='headaches')
      # Generate GAMs using the formula that incorporates scan reason
      gamm <- createMainGamm(df, phenotype)
  
      # Predict on the GAMs for the actual data
      gammPreds <- predict_gam(gamm$gam, values = modelFixedValues)
    }

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
regionalPeaksM <- generateHemispherePlots(greyVolDf, regionalPhenotypes, parsedRegionalPhenotypes, '', "M")
regionalPeaksF <- generateHemispherePlots(greyVolDf, regionalPhenotypes, parsedRegionalPhenotypes, '', "F")

## Step 11: Compare CLIP age at peak to Lifespan age at peak -------------------
# Read the ages from lifespan
lifespanDf <- read.csv('/Users/youngjm/Data/lifespan_growth_charts/Lifespan_Data_Peaks_Table_2_2.csv')
lifespanDf <- merge(regionalPeaks, lifespanDf, by='feat')
regionalPeaks$Peak <- as.numeric(regionalPeaks$em)
regionalPeaksM$Peak <- as.numeric(regionalPeaksM$em)
regionalPeaksF$Peak <- as.numeric(regionalPeaksF$em)

lifespanDf$peakDiff <- as.numeric(lifespanDf$Peak) - as.numeric(regionalPeaks$Peak)
regionMin <- min(regionalPeaks$Peak, regionalPeaksM$Peak, regionalPeaksF$Peak, lifespanDf$Peak)
regionMax <- max(regionalPeaks$Peak, regionalPeaksM$Peak, regionalPeaksF$Peak, lifespanDf$Peak)

regionalPeaksPlot <- regionalPeaks %>% 
  ggseg(mapping=aes(fill=Peak),
        hemisphere='left') 

regionalPeaksMPlot <- regionalPeaksM %>% 
  ggseg(mapping=aes(fill=Peak),
        hemisphere='left') 

regionalPeaksFPlot <- regionalPeaksF %>% 
  ggseg(mapping=aes(fill=Peak),
        hemisphere='left') 

lifespanPeaksPlot <- lifespanDf %>% 
  ggseg(mapping=aes(fill=Peak),
        hemisphere='left') 

diffBarPlot <- ggplot(data=lifespanDf, aes(y=region, x=abs(peakDiff), fill=region)) +
  geom_bar(stat = 'identity', position = 'identity') +
  # coord_flip() +
  theme_minimal() +
  guides(fill = 'none') +
  ylab('Region') + xlab('Difference in Age at Peak from Lifespan')


# Plot 3 graphs in 1 figure
regionalPlots <- (regionalPeaksPlot + lifespanPeaksPlot) /
  (regionalPeaksMPlot + regionalPeaksFPlot) +
  plot_layout(guides = "collect") & 
  scale_fill_viridis_c(limits=c(0,regionMax))
regionalPlots

# (test) + diffBarPlot + plot_layout(layout)

## Calculate the correlation between Lifespan and CLIP age at peak -------------
tmp <- lifespanDf[regionalPeaks$feat != 'medialorbitofrontal',]

# for each phenotype, get the average volume from greyVolDf
avgVol <- c()
for (pheno in tmp$feat){
  print(pheno)
  avg <- mean(greyVolDf[[pheno]])/4000
  avgVol <- append(avgVol, avg)
}
tmp$avgVol <- avgVol

minVal <- min(tmp$Peak, as.numeric(tmp$em))
maxVal <- max(tmp$Peak, as.numeric(tmp$em))

ggplot(data=tmp, aes(color=as.factor(feat), shape=as.factor(feat), fill=as.factor(feat)))+
  geom_point(aes(x=as.numeric(em), y=Peak), size=avgVol, alpha=0.65) +
  geom_abline(slope = 1) +
  scale_shape_manual(values = rep(21:25, 6)) +
  scale_color_manual(values = rep(cbbPalette, 4)) +
  scale_fill_manual(values = rep(cbbPalette, 4)) +
  # theme_minimal() +
  guides(guide_legend(title="Region")) +
  # expand_limits(x=c(0, max(regionalPeaks$em)), y=c(0, max(regionalPeaks$Peak))) +
  xlab('Age at Peak CLIP (years)') +
  ylab('Age at Peak Lifespan (years)') +
  xlim(minVal-0.5, maxVal+0.5) +
  ylim(minVal-0.5, maxVal+0.5) +
  labs(title = 'Age at Peak: Lifespan vs. CLIP') 

cor(tmp$Peak, as.numeric(tmp$em))

## Step 6: Load SynthSeg outputs -----------------------------------------------

synthsegFn <- '/Users/youngjm/Data/clip/images/derivatives/metrics_synthsegplus.csv'
synthsegDf <- read.csv(synthsegFn)

# Drop any rows where neurofibromatosis is in the scan_reason_categories column
synthsegDf <- synthsegDf[!grepl("neurofibromatosis", synthsegDf$scan_reason_primary), ]
synthsegDf <- synthsegDf[!grepl("neurofibromatosis", synthsegDf$scan_reason_categories), ]
synthsegDf <- synthsegDf %>%
  mutate(confirm_neurofibromatosis = case_when(
    grepl('1', confirm_neurofibromatosis, fixed=TRUE) ~ TRUE,
    grepl('0', confirm_neurofibromatosis, fixed=TRUE) ~ FALSE,
    grepl('False', confirm_neurofibromatosis, fixed=TRUE) ~ FALSE,
    grepl('', confirm_neurofibromatosis, fixed=TRUE) ~ FALSE,
    is.na(confirm_neurofibromatosis) ~ FALSE
  ))
synthsegDf <- synthsegDf[(synthsegDf$confirm_neurofibromatosis == FALSE), ]

# Make an age in years column
synthsegDf$age_in_years <- synthsegDf$age_in_days/365.25

synthsegDf$Processing <- "SynthSeg+"

# Drop any 1.5T scans
synthsegDf <- synthsegDf[synthsegDf$MagneticFieldStrength != "1.5",]

# Drop any scans with ratings less than 0
synthsegDf <- synthsegDf[synthsegDf$rawdata_image_grade >= 0, ]

# Rename columns - WRONG
names(synthsegDf)[names(synthsegDf) == "pat_id"] <- "patient_id"

# Only one scan per subject
# Sort the dataframe by patient_id and scanner_id
synthsegDf <- synthsegDf[ with(synthsegDf, order(synthsegDf$patient_id, synthsegDf$scan_id)), ]
synthsegDf <- synthsegDf[!duplicated(synthsegDf$patient_id), ]
synthsegDf <- synthsegDf[complete.cases(synthsegDf), ]
synthsegDf$patient_id <- droplevels(as.factor(synthsegDf$patient_id))

# Add new column: going to sum TotalGrayVol + CerebralWhiteMatterVol
synthsegDf$TCV <- synthsegDf$SS_TotalGrayVolume + synthsegDf$SS_CerebralWhiteMatterVolume

# Some of the columns in this df should be factors
toFactor <- c('sex', 'Processing', 'MagneticFieldStrength', 'scanner_id',
              'scan_reason_primary')
synthsegDf[toFactor] <- lapply(synthsegDf[toFactor], factor)
synthsegDf <- synthsegDf[(synthsegDf$rawdata_image_grade >= 1),]
synthsegDf <- addPrimaryScanReasonCol(synthsegDf)

tableColumns <- c("Total Gray\nVolume", "Cerebral White\nMatter Volume", 
                  "Subcortical Gray\nVolume", "Ventricle\nVolume", 
                  "Total Cerebrum\nVolume") 

tableRows <- c("Surface Holes", 
               "Sex = Male", 
               "Age")

miniGlobalCols <- c("SS_TotalGrayVolume", "SS_CerebralWhiteMatterVolume", "SS_SubcorticalGrayVolume",
                    "SS_VentricleVolume", "TCV")

synthsegPlotsTables <- generatePlotsAndTablesForDataset(synthsegDf, miniGlobalCols, tableColumns, tableRows, 'SynthSeg+ Phenotypes')





## Step ???: PCA?

# First regression to remove the effects of age and sex
regressOutAgeSex <- function(phenotype, df){
  formula <- as.formula(paste(phenotype, "age_in_years + sex + SurfaceHoles", sep="~"))
  model <- gam(formula, data=df)
  return(model$residuals)
}

# Run regression on the global phenotypes
residualGlobalPhenoCols <- c()
for (pheno in globalPhenoCols){
  newPheno <- paste0(pheno, "_residuals")
  highQDf[[newPheno]] <- regressOutAgeSex(pheno, highQDf)
  combattedHighQDf[[newPheno]] <- regressOutAgeSex(pheno, combattedHighQDf)
  residualGlobalPhenoCols <- append(residualGlobalPhenoCols, newPheno)
}

# Do the PCA
pcCols <- c(residualGlobalPhenoCols)
pc <- prcomp(highQDf[, pcCols], center=TRUE, scale.=TRUE)
pcCombatted <- prcomp(combattedHighQDf[, pcCols], center=TRUE, scale.=TRUE)
pcCentiles <- prcomp()

# PCA non-combatted
summary(pc)

plotPCsForGroup <- function(pc, group){
  a12 <- ggbiplot(pc,
                  choices = 1:2,
                  alpha = 0.35,
                  ellipse = TRUE,
                  groups = group)
  a13 <- ggbiplot(pc,
                  choices = c(1,3),
                  alpha = 0.35,
                  ellipse = TRUE,
                  groups = group)
  a14 <- ggbiplot(pc,
                  choices = c(1,4),
                  alpha = 0.35,
                  ellipse = TRUE,
                  groups = group)
  a23 <- ggbiplot(pc,
                  choices = c(2,3),
                  alpha = 0.35,
                  ellipse = TRUE,
                  groups = group)
  a24 <- ggbiplot(pc,
                  choices = c(2,4),
                  alpha = 0.35,
                  ellipse = TRUE,
                  groups = group)
  a34 <- ggbiplot(pc,
                  choices = c(3,4),
                  alpha = 0.35,
                  ellipse = TRUE,
                  groups = group)
  (a12 | a13 | a14) / (a23 | a24 | a34) + plot_layout(guides = "collect")
}

# Plots for non-combatted residuals
plotPCsForGroup(pc, highQDf$top_scan_reason_factors)
plotPCsForGroup(pc, highQDf$scanner_id)
c <- ggscreeplot(pc)
c

# Plots for combatted residuals
plotPCsForGroup(pcCombatted, combattedHighQDf$top_scan_reason_factors)
plotPCsForGroup(pcCombatted, combattedHighQDf$scanner_id)


# d <- dist(highQDf[, pcCols]) # euclidean distances between the rows
# fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# fit # view results
# # plot solution
# x <- fit$points[,1]
# y <- fit$points[,2]
# 
# ggplot(data=highQDf, aes(x=x, y=y, fill=.data$scanner_id)) + geom_density2d()
# geom_density2d(aes(x=x, y=y))

# PCA of the residuals from the combatted data
summary(pcCombatted)
a <- ggbiplot(pcCombatted,
              alpha = 0.5,
              groups = combattedHighQDf$top_scan_reason_factors)

b <- ggbiplot(pcCombatted,
              alpha = 0.5,
              groups = combattedHighQDf$scanner_id) 

c <- ggscreeplot(pcCombatted)

(a / b) | c
## TESTING ---------------------------------------------------------------------
# Load the expected centiles from Jakob
predCentFn <- "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/2022-06-02_ages_sex_for_lifespand_prediction_GMV.csv"
predGMVDf <- read.csv(predCentFn)
predGMVDf$sex <- as.factor(predGMVDf$sex)

gmvCor <- cor(combattedDf$TotalGrayVol, predGMVDf$life_GMV)
gmvCor

# Scatterplot: CLIP vs. lifespan-predicted GMV
scatteredCombat <- ggplot() +
  geom_point(aes(x=combattedDf$TotalGrayVol, y=predGMVDf$life_GMV), alpha=0.65) +
  geom_abline(slope = 1) +
  # geom_smooth(method='lm', formula= predGMVDf$life_GMV~combattedDf$TotalGrayVol, color='red') +
  # scale_shape_manual(values = rep(21:25, 6)) +
  # scale_color_manual(values = rep(cbbPalette, 4)) +
  # scale_fill_manual(values = rep(cbbPalette, 4)) +
  # theme_minimal() +
  # guides(guide_legend(title="Region")) +
  # expand_limits(x=c(0, max(regionalPeaks$em)), y=c(0, max(regionalPeaks$Peak))) +
  xlab('PostCombat GMV for CLIP') +
  ylab('GMV Predicted from Lifespan') +
  # xlim(minVal-0.5, maxVal+0.5) +
  # ylim(minVal-0.5, maxVal+0.5) +
  labs(title = 'Lifespan-Predicted GMV vs. PostCombat CLIP GMV') 

scatteredCC <- ggplot() +
  geom_point(aes(x=ccDf$GMVTransformed.normalised, y=predGMVDf$life_GMV), alpha=0.65) +
  geom_abline(slope = 10000) +
  xlab('Normalised PostCombat GMV for CLIP') +
  ylab('GMV Predicted from Lifespan') +
  labs(title = 'Lifespan-Predicted GMV vs. Normalised PostCombat CLIP GMV') 

# a <- combatPlotsTables$predictionPlots$TotalGrayVol + 
#   xlab("Age (years)") +
#   ylab("Total Gray Volume (mm3)") +
#   title("PostCombat")
#   
# b <- centiledCombattedPlotsTables$predictionPlots$GMVTransformed.normalised +
#   xlab("Age (years)") +
#   ylab("Total Gray Volume (mm3)") +
#   title("Normalised PostCombat")

(combatPlotsTables$predictionPlots$TotalGrayVol / centiledCombattedPlotsTables$predictionPlots$GMVTransformed.normalised) | (scatteredCombat / scatteredCC )
cor.test(x=combattedDf$TotalGrayVol, y=predGMVDf$life_GMV)
cor.test(x=ccDf$GMVTransformed.normalised, y=predGMVDf$life_GMV)


## TESTING: Regional Trajectories ----------------------------------------------




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
  gammScanReason <- createMainGamm(analysisDf, phenotype)
  
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