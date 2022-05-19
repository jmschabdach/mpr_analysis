gc()
dev.off(dev.list()["RStudioGD"])
# Libraries --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgcv)
library(tidymv)
library(patchwork) # graph organization within a figure
library(gtsummary)
library(grid)
library(stringr)
library(gridExtra)
library(reshape2)
library(tables)
library(grid)
library(gridExtra)
library(ggseg)
library(data.table)
library(formattable)
library(tidyr)

# Colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Diverging palette with 6 classes
assign("myPalette", c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"),
       envir = .GlobalEnv)
assign("fig_pub_theme", theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_line(color = 'gray'), #remove major gridlines
  # panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
), envir = .GlobalEnv)

# Specify the file names
clipFsFn <- '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_7.1.1_tables/mpr_fs_reconall_7.1.1_structural_stats.csv'
clipIfsFn <- '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_7.1.1_tables/mpr_ifs_reconall_7.1.1_structural_stats.csv'
# clipFn <- '/Users/youngjm/Data/clip/tables/CLIPv0.7/2022-03_analysis_features.csv'
q22FsFn <- '/Users/youngjm/Data/22q11/tables/2022-04_fs_demo_phenotypes.csv'
q22IfsFn <- '/Users/youngjm/Data/22q11/tables/2022-04_ifs_demo_phenotypes.csv'
enigmaFn <- "/Users/youngjm/Data/enigma_22q/2022-04-13_combined_demographics_phenotypes.csv"

# FUNCTIONS --------------------------------------------------------------------
prepInfo <- function(data.df){
  # Adding column based on other column:
  new.data.df <- data.df %>%
    mutate(MagneticFieldStrength = case_when(
      grepl('1p5', scan_id, fixed=TRUE) ~ "1.5",
      grepl('3p0', scan_id, fixed=TRUE) ~ "3.0"
    ))
  # Make the new column a factor
  new.data.df$MagneticFieldStrength <- as.factor(new.data.df$MagneticFieldStrength)
  # Wait let's also make the scanner_id column a factor here too
  new.data.df$scanner_id <- as.factor(new.data.df$scanner_id)
  # And add an age in years column
  new.data.df$age_in_years <- new.data.df$age_at_scan_days/365.25
  # Let's also drop all 1.5 scans
  new.data.df <- new.data.df[new.data.df$MagneticFieldStrength == "3.0", ]
  
  return(new.data.df)
}

createGamm <- function(df, measure) {
  formula <- as.formula(paste(measure, "s(log(age_in_years), fx=T) +
                      # s(log(age_in_years), by = ordered(Diagnosis), fx=T) +
                      SurfaceHoles +
                      ordered(Diagnosis) +
                      sex", sep="~"))

  mixedModel <- gamm(formula,
                     random = list(scanner_id=~1),
                     data = df)
  # Return the model
  return(mixedModel)
}

createDiagnosisGamm <- function(df, measure) {
  formula <- as.formula(paste(measure, "s(log(age_in_years), fx=T) +
                      SurfaceHoles +
                      sex", sep="~"))
  
  mixedModel <- gamm(formula,
                     random = list(scanner_id=~1),
                     data = df)
  # Return the model
  return(mixedModel)
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
# Make a scatter plot of a specific measure for the original data. Include the
# confidence intervals that come from the output of the GAMM predictions
# @param predictions Output from predict_gam() for the measure and its gam
# @param origData A dataframe containing the original data
# @param measure A column of the original dataframe
# @param measureTitle A string used to describe the measure in the plot
# @returns plot01 A ggplot object with scatter points and a confidence interval
generatePlotScatter <- function(predictions, origData, measure, measureTitle, dx){
  
  peakAge <- getAgeAtPeak(predictions)
  
  # Filter by diagnosis
  if (dx == "HC"){
    # Colorblind palette with black:
    cbbPalette <- myPalette[6:4]
    origData <- origData[origData$Diagnosis == 'HC',]
  } else {
    # Colorblind palette with black:
    cbbPalette <- myPalette[1:3]
    origData <- origData[origData$Diagnosis == '22q11DS',]
  }
  
  if (grepl('V', measureTitle)) {
    ylabel=expression(paste("Volume (", 'mm'^3, ')'))
  } else if (grepl('Area', measureTitle)) {
    ylabel=expression(paste("Surface Area (", 'mm'^2, ')'))
  } else if (grepl('Thickness', measureTitle)) {
    ylabel='Thickness (mm)'
  } else {
    ylabel=''
  }
  
  measure <- origData[ , measure]
  predictions <- predictions[predictions$Diagnosis == dx, ]
  
  plot01 <- predictions %>%
    ggplot(aes(x=age_in_years, y=fit)) +
    geom_point(data=origData, aes(x=age_in_years, y=measure, color=Group), alpha=0.5) +
    geom_smooth_ci() +
    scale_color_manual(values = cbbPalette, name = "Diagnosis and Processing") +
    # theme(axis.title = element_blank()) +
    fig_pub_theme +
    labs(title = measureTitle,
         x = 'Age (years)',
         y = ylabel)
  # print(plot01)
  
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
  # rownames(subTab) <- c('Age', 'Age by Dx')
  
  
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
                   'Diagnosis',
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
# Make a table showing the statistical significance of linear parameters on image phenotype
# @param aGam A gam object
# @param nTests A number indicating how many tests are being performed/should be corrected for
# @returns linT A table object to be displayed using gridExtra
getBetasFromGamSummary <- function(aGam){
  
  # Get the p-values
  pvals <- signif(summary(aGam)$p.pv, 3)
  
  serrs <- signif(summary(aGam)$se, 3)
  tvals <- signif(summary(aGam)$p.t, 3)
  betas <- c()
  
  for (i in names(pvals)){
    betas[[i]] <- (tvals[[i]] * serrs[[i]])
  }
  
  pvals[['Age']] <- signif(summary(aGam)$s.pv, 3)
  # print(pvals)
  subTab <- as.data.frame(betas)

  return(subTab)
}

getDiagnosisStats <- function(aGam){
  subTab <- signif(summary(aGam)$p.table, 3)
  
  newRowNames <- c('Surface Holes',
                   'Diagnosis',
                   'Sex')
  rows <- rownames(subTab)[2:length(rownames(subTab))]
  subTab <- subset(subTab, rownames(subTab) %in% rows)
  rownames(subTab) <- newRowNames
  
  # Rename and trim the columns
  colnames(subTab) <- c('beta', 'std err', 't-value', 'p-value')
  cols <- c('p-value', 'beta')
  subTab <- subTab['Diagnosis',cols]
}

generateDiagnosisPlotCI <- function(predictions, origData, measure, measureTitle){
  
  dxColors <- c(myPalette[1], myPalette[6])
  
  if (grepl('V', measureTitle)) {
    ylabel=expression(paste("Volume (", 'mm'^3, ')'))
  } else if (grepl('Area', measureTitle)) {
    ylabel=expression(paste("Surface Area (", 'mm'^2, ')'))
  } else if (grepl('Thickness', measureTitle)) {
    ylabel='Thickness (mm)'
  } else {
    ylabel=''
  }
  
  plot02 <- predictions %>%
    ggplot(aes(x=age_in_years, y=fit, color=Diagnosis)) +
    geom_smooth_ci() +
    scale_color_manual(values = dxColors) +
    # scale_x_continuous(trans='log10') +
    # theme(axis.title = element_blank()) +
    fig_pub_theme +
    labs(title = measureTitle,
         x = 'Age (years)',
         y = ylabel)
  
  return(plot02)
}

makeTablePretty <- function(secondStep, tableColumns, title){
  green0 <- "#4fa64f"
  green1 <- "#85c285" 
  green2 <- "#abd5ab"
  significance_formatter <- function (...) {
    formatter("span", style = function(x) {
      style(display = "block",
            padding = "0 4px", 
            `border-radius` = "4px", 
            `background-color` = ifelse(x < 0.0005 , green0, 
                                        ifelse(x < 0.005, green1, 
                                               ifelse(x < 0.05, green2, 'white')
                                        )
            )
      ) # Remember to change the colors!
    })}
  
  formatted <- list()
  for (col in tableColumns){
    formatted[[col]] <- significance_formatter()
  }
  
  obj <- formattable(secondStep,
              align =c("c","c","c","c", "c", "c", "c"), 
              formatted,
              caption=title)
  obj
  return(obj)
}

generateHorizontalTrajectoriesAndTables <- function(myData, phenotypes, tableColumns, tableRows, title){
  # Initialize variables
  # Set up empty variables to generate table later
  ageAtPeak <- list()
  predictionPlots <- list()
  tableValues <- list()
  tableBetas <- list()
  
  modelFixedValues <- list(SurfaceHoles = mean(myData$SurfaceHoles), 
                           sex='M')
  
  hcTitles <- c("Total Gray Volume (mm^3)", "Cerebral White Matter Volume (mm^3)",
                "Subcortical Gray Matter Volume (mm^3)", "eTIV (mm^3)", "Ventricle Volume (mm^3)",
                "Cortical Surface Area (mm^2)", "Mean Cortical Thickness (mm)")
  otherTitles <- c("", "", "", "", "", "", "")
  # for phenotype in phenotypes...
  for (i in 1:length(phenotypes)){
    phenotype <- phenotypes[[i]]
    print(phenotype)
    # Generate GAMs using the formula that incorporates scan reason
    gammScanReason <- createGamm(myData, phenotype)
    
    # Predict on the GAMs for the actual data
    gammScanReasonPreds <- predict_gam(gammScanReason$gam, values = modelFixedValues)
    
    # Get the age at peak
    ageAtPeak[[phenotype]] <- getAgeAtPeak(gammScanReasonPreds)
    
    # Generate scatter plots with confidence intervals
    plotScatterScanReason <- generatePlotScatter(gammScanReasonPreds, myData, phenotype, hcTitles[[i]], 'HC')
    predictionPlots[[paste(phenotype, 'scatter_hc', sep='_')]] <- plotScatterScanReason
    plotScatterScanReason <- generatePlotScatter(gammScanReasonPreds, myData, phenotype, otherTitles[[i]], '22q11DS')
    predictionPlots[[paste(phenotype, 'scatter_22q', sep='_')]] <- plotScatterScanReason
    predictionPlots[[paste(phenotype, 'ci', sep='_')]] <- generateDiagnosisPlotCI(gammScanReasonPreds, myData, phenotype, otherTitles[[i]])
    
    tParam <- generateScanReasonsParametricTable(gammScanReason$gam)
    tLin <- generateScanReasonsLinearTable(gammScanReason$gam)
    tableValues[[phenotype]] <- append(tLin, tParam)
    tableBetas[[phenotype]] <- getBetasFromGamSummary(gammScanReason$gam)
  }
  
  # Plot all scatter/CI plots in 1 figure
  hcRow <- names(predictionPlots)[grepl('hc', names(predictionPlots))]
  q22Row <- names(predictionPlots)[grepl('22q', names(predictionPlots))]
  ciRow <- names(predictionPlots)[grepl('ci', names(predictionPlots))]
  hcPatchwork <- wrap_plots(predictionPlots[hcRow], nrow=1, guides="collect") 
  q22Patchwork <- wrap_plots(predictionPlots[q22Row], nrow=1, guides="collect") 
  ciPatchwork <- wrap_plots(predictionPlots[ciRow], nrow=1, guides="collect") 
  
  # 2800 x 700
  p <- wrap_plots(predictionPlots, nrow=3, guides = "collect", byrow=FALSE)
  print(hcPatchwork / q22Patchwork / ciPatchwork)
  
  # Make table: p-values of factor/phenotype
  firstStep <- lapply(tableValues, unlist)
  secondStep <- as.data.frame(firstStep, stringsAsFactors = F)
  colnames(secondStep) <- tableColumns
  rownames(secondStep) <- tableRows
  t1 <- makeTablePretty(secondStep, tableColumns, paste0("Significance of Features for", title, " (p-values)", sep=' '))
  print(t1)
  
  # Make table: betas
  firstStep <- lapply(tableBetas, unlist)
  secondStep <- as.data.frame(firstStep, stringsAsFactors = F)
  colnames(secondStep) <- tableColumns
  # rownames(secondStep) <- tableRows[1:(length(tableRows)-1)]
  t2 <- makeTablePretty(secondStep, tableColumns, paste0("Significance of Features for ", title, " (beta values)", sep=' '))
  print(t2)
}

generateVerticalTrajectories <- function(myData, phenotypes, tableColumns, tableRows, title){
  # Initialize variables
  # Set up empty variables to generate table later
  ageAtPeak <- list()
  predictionPlots <- list()
  tableValues <- list()
  tableBetas <- list()
  
  modelFixedValues <- list(SurfaceHoles = mean(myData$SurfaceHoles), 
                           sex='M')
  
  hcTitles <- c("Total Gray Volume", "Cerebral White Matter Volume",
                "Subcortical Gray Matter Volume", "eTIV", "Ventricle Volume",
                "Cortical Surface Area", "Mean Cortical Thickness")
  otherTitles <- c("", "", "", "", "", "", "")
  # for phenotype in phenotypes...
  for (i in 1:length(phenotypes)){
    phenotype <- phenotypes[[i]]
    print(phenotype)
    # Generate GAMs using the formula that incorporates scan reason
    gammScanReason <- createGamm(myData, phenotype)
    
    # Predict on the GAMs for the actual data
    gammScanReasonPreds <- predict_gam(gammScanReason$gam, values = modelFixedValues)
    
    # Get the age at peak
    ageAtPeak[[phenotype]] <- getAgeAtPeak(gammScanReasonPreds)
    
    # Generate scatter plots with confidence intervals
    plotScatterScanReason <- generatePlotScatter(gammScanReasonPreds, myData, phenotype, hcTitles[[i]], 'HC')
    predictionPlots[[paste(phenotype, 'scatter_hc', sep='_')]] <- plotScatterScanReason
    plotScatterScanReason <- generatePlotScatter(gammScanReasonPreds, myData, phenotype, otherTitles[[i]], '22q11DS')
    predictionPlots[[paste(phenotype, 'scatter_22q', sep='_')]] <- plotScatterScanReason
    predictionPlots[[paste(phenotype, 'ci', sep='_')]] <- generateDiagnosisPlotCI(gammScanReasonPreds, myData, phenotype, otherTitles[[i]])
    
    tParam <- generateScanReasonsParametricTable(gammScanReason$gam)
    tLin <- generateScanReasonsLinearTable(gammScanReason$gam)
    tableValues[[phenotype]] <- append(tLin, tParam)
    tableBetas[[phenotype]] <- getBetasFromGamSummary(gammScanReason$gam)
  }
  
  # Plot all scatter/CI plots in 1 figure
  # hcRow <- names(predictionPlots)[grepl('hc', names(predictionPlots))]
  # q22Row <- names(predictionPlots)[grepl('22q', names(predictionPlots))]
  # ciRow <- names(predictionPlots)[grepl('ci', names(predictionPlots))]
  # hcPatchwork <- wrap_plots(predictionPlots[hcRow], ncol=1, guides="collect") 
  # q22Patchwork <- wrap_plots(predictionPlots[q22Row], ncol=1, guides="collect") 
  # ciPatchwork <- wrap_plots(predictionPlots[ciRow], ncol=1, guides="collect") 
  
  # 2800 x 700
  p <- wrap_plots(predictionPlots, ncol=3, guides = "collect", byrow=TRUE)
  # print(hcPatchwork / q22Patchwork / ciPatchwork)
  print(p)
}

generateTableScaledBetas <- function(myData, phenotypes, tableColumns, tableRows, title){
  tableValues <- list()
  tableBetas <- list()
  modelFixedValues <- list(SurfaceHoles = mean(myData$SurfaceHoles), 
                           sex='M')
  # for phenotype in phenotypes...
  for (i in 1:length(phenotypes)){
    phenotype <- phenotypes[[i]]
    print(phenotype)
    # Generate GAMs using the formula that incorporates scan reason
    myData[[phenotype]] <- scale(myData[[phenotype]])
    gammScanReason <- createGamm(myData, phenotype)
    
    # # Predict on the GAMs for the actual data
    # gammScanReasonPreds <- predict_gam(gammScanReason$gam, values = modelFixedValues)
    
    # dxStats[[phenotype]] <- getDiagnosisStats(gammScanReason$gam)
    tParam <- generateScanReasonsParametricTable(gammScanReason$gam)
    tLin <- generateScanReasonsLinearTable(gammScanReason$gam)
    print(tParam)
    print(tLin)
    tableValues[[phenotype]] <- append(tLin, tParam)
    tableBetas[[phenotype]] <- getBetasFromGamSummary(gammScanReason$gam)
  }
  
  # Make table: diagnosis stats
  firstStep <- lapply(tableValues, unlist)
  print(firstStep)
  secondStep <- as.data.frame(firstStep, stringsAsFactors = F)
  colnames(secondStep) <- tableColumns
  t3 <- makeTablePretty(secondStep, tableColumns, paste0("Scaled Significance, p-values, ,", title, sep=' '))
  # print(t3)
  
  # Make table: diagnosis stats
  firstStep <- lapply(tableBetas, unlist)
  secondStep <- as.data.frame(firstStep, stringsAsFactors = F)
  colnames(secondStep) <- tableColumns
  t4 <- makeTablePretty(secondStep, tableColumns, paste0("Scaled Significance, betas,", title, sep=' '))
  t4
  
}

generateAgeSexDistributionPlot <- function(dataDf, diagnosis, colorPalette, title, maxAge){
  # Set the title string
  if (diagnosis == 'HC'){
    title = paste(title, "Control Age at Scan", sep=' ')
    dataDf <- dataDf[dataDf$Diagnosis == 'HC', ]
    cols <- colorPalette[6:5]
  } else {
    title = paste(title, "22q11DS Age at Scan", sep=' ')
    dataDf <- dataDf[dataDf$Diagnosis == '22q11DS', ]
    cols <- colorPalette[1:2]
  }
  print(diagnosis)
  newPlot <- ggplot(data=dataDf, aes(x=age_in_years, fill=as.factor(sex))) +
    geom_histogram(position="stack", binwidth=1) +
    scale_fill_manual(values = cols, drop=FALSE) + 
    xlim(0, maxAge) +
    labs(title = title,
         x = 'Subject Age (Years)',
         y = '# Subjects',
         fill="Sex")
  
  return(newPlot)
}

getRegionalPhenotypes <- function(brainDf, metadataCols){
  # Get the columns of the dataframe
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
  
  # Get all columns in the dataframe that end with "_grayVol"
  greyVolCols <- brainDfCols[endsWith(brainDfCols, '_grayVol')]
  
  # Get all columns in the list of grayVol columns with partial string matches to the list of regional phenotypes
  phenoCols <- c()
  for (pheno in parsedRegionalPhenotypes) {
    cols <- greyVolCols[grepl(paste("_",pheno, sep=''), greyVolCols, fixed=TRUE)]
    phenoCols <- c(phenoCols, cols)
  }
  
  # Separate the left and right hemisphere columns
  grayVolLhCols <- sort(phenoCols[startsWith(phenoCols, 'lh_')])
  grayVolRhCols <- sort(phenoCols[startsWith(phenoCols, 'rh_')])
  
  # Let's calculate the bilateral average of each regional phenotype
  for (i in 1:length(parsedRegionalPhenotypes)){
    brainDf[[parsedRegionalPhenotypes[[i]]]] <- (brainDf[[grayVolLhCols[[i]]]] + brainDf[[grayVolRhCols[[i]]]])/2
  }
  
  greyVolDf <- brainDf[append(metadataCols, parsedRegionalPhenotypes)]
  return(greyVolDf)
}

##
# Make the composite plot and tables for a given dataframe/set of phenotypes
# @param df A dataframe where rows are scans and columns are phenotypes and features
# @param phenotypes A vector of phenotypes to analyze
# @param colNames A vector of cleaned names to use for the table columns
# @param rowNames A vector of cleaned names to use for the table rows
# @param title A string to differentiate the figures for this dataframe vs. other dataframes
generateHemispherePlots <- function(df, title){
  # Initialize variables
  phenotypes <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal',
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
  phenotypesNoWS <- c('bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal',
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
  # Set up empty variables to generate table later
  ageAtPeak <- c()
  
  modelFixedValues <- list(SurfaceHoles = mean(myData$SurfaceHoles), 
                           sex='M')
  
  # for phenotype in phenotypes...
  for (phenotype in phenotypesNoWS){
    # Generate GAMs using the formula that incorporates scan reason
    gamm <- createDiagnosisGamm(df, phenotype)
    
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
  
  # print(results)
  
  # Plot the dataframe using ggseg
  p <- results %>% 
    ggseg(mapping=aes(fill=as.numeric(em)),
          hemisphere='left') +
    labs(title = paste(title, "Age at Peak (years)"), legend="Age") +
    theme(axis.title = element_blank()) +
    scale_fill_gradient(low = "blue", high = "red", na.value = NA) 
  
  grid.arrange(p)
  return(list(results, p))
}

addPrimaryScanReasonCol <- function(df){
  if ("scan_reason_primary" %in% colnames(df)){
    # get the top 6 primary reasons
    topScanReasons <- names(sort(table(df$scan_reason_primary), decreasing=TRUE)[1:5])
    
    # build a column with these top 5 primary reasons and a generous other category
    df <- mutate(df, top_scan_reason_factors = if_else(is.element(scan_reason_primary, topScanReasons),
                                                       paste(scan_reason_primary),
                                                       "other"))
    df$top_scan_reason_factors <- as.factor(df$top_scan_reason_factors)
    # Put the "other" category first
    df$top_scan_reason_factors <- relevel(df$top_scan_reason_factors, "other")
  } else {
    df$top_scan_reason_factors <- NA
  }
  
  return(df)
}

# Step 1: load and prep data ---------------------------------------------------

# Read the data from the files
clipFsDf <- read.csv(clipFsFn, stringsAsFactors = TRUE)
clipIfsDf <- read.csv(clipIfsFn, stringsAsFactors = TRUE)
clin22qFsDf <- read.csv(q22FsFn, stringsAsFactors = TRUE)
clin22qIfsDf <- read.csv(q22IfsFn, stringsAsFactors = TRUE)

# Drop any clip NF1
clipFsDf <- clipFsDf[!grepl("neurofibromatosis", clipFsDf$scan_reason_primary), ]
clipFsDf <- clipFsDf[!grepl("neurofibromatosis", clipFsDf$scan_reason_categories), ]
clipFsDf <- clipFsDf[!grepl("true", clipFsDf$confirm_neurofibromatosis), ]
clipIfsDf <- clipIfsDf[!grepl("neurofibromatosis", clipIfsDf$scan_reason_primary), ]
clipIfsDf <- clipIfsDf[!grepl("neurofibromatosis", clipIfsDf$scan_reason_categories), ]
clipIfsDf <- clipIfsDf[!grepl("true", clipIfsDf$confirm_neurofibromatosis), ]

# Making sure column names are the same
names(clin22qFsDf)[names(clin22qFsDf) == 'scanId'] <- 'scan_id'
names(clin22qIfsDf)[names(clin22qIfsDf) == 'scanId'] <- 'scan_id'
names(clin22qFsDf)[names(clin22qFsDf) == 'subjId'] <- 'patient_id'
names(clin22qIfsDf)[names(clin22qIfsDf) == 'subjId'] <- 'patient_id'

clipFsDf <- prepInfo(clipFsDf)
clipIfsDf <- prepInfo(clipIfsDf)
clin22qFsDf <- prepInfo(clin22qFsDf)
clin22qIfsDf <- prepInfo(clin22qIfsDf)

# Add a Diagnosis column to each data frame to specify the group the data belongs to
clipFsDf$Diagnosis <- 'HC'
clipIfsDf$Diagnosis <- 'HC'
clin22qFsDf$Diagnosis <- '22q11DS'
clin22qIfsDf$Diagnosis <- '22q11DS'

# Add a Processing column to each data frame to specify the group the data belongs to
clipFsDf$FSVersion <- 'FS71'
clipIfsDf$FSVersion <- 'IFS'
clin22qFsDf$FSVersion <- 'FS71'
clin22qIfsDf$FSVersion <- 'IFS'

# Remove bad rawdata grades
clipFsDf <- clipFsDf[clipFsDf$rawdata_image_grade > 0, ]
clipIfsDf <- clipIfsDf[clipIfsDf$rawdata_image_grade > 0, ]

# Load ENIGMA data
enigmaDf <- read.csv(enigmaFn)

# Make a age in years column
enigmaDf$age_in_years <- enigmaDf$age_days/365.25

# Drop 22q11Dup
enigmaDf <- enigmaDf[!grepl("22q11Dup", enigmaDf$dx), ]
enigmaDf <- enigmaDf[, -which(grepl("rh", colnames(enigmaDf)))]
enigmaDf <- enigmaDf[, -which(grepl("lh", colnames(enigmaDf)))]
enigmaDf <- enigmaDf[, -which(grepl("Right", colnames(enigmaDf)))]
enigmaDf <- enigmaDf[, -which(grepl("Left", colnames(enigmaDf)))]

# Rename
names(enigmaDf)[names(enigmaDf) == 'subjId'] <- 'patient_id'
names(enigmaDf)[names(enigmaDf) == 'scanId'] <- 'scan_id'
names(enigmaDf)[names(enigmaDf) == 'fs_version'] <- 'FSVersion'
names(enigmaDf)[names(enigmaDf) == 'dx'] <- 'Diagnosis'
names(enigmaDf)[names(enigmaDf) == 'Ventricles'] <- 'VentricleVolume'
names(enigmaDf)[names(enigmaDf) == 'WMV'] <- 'CerebralWhiteMatterVol'
names(enigmaDf)[names(enigmaDf) == 'site'] <- 'scanner_id'
names(enigmaDf)[names(enigmaDf) == 'age_days'] <- 'age_at_scan_days'

# Change the enigma sexes to M/F
enigmaDf <- enigmaDf %>%
  mutate(sex = case_when(
    grepl('Male', sex, fixed=TRUE) ~ "M",
    grepl('Female', sex, fixed=TRUE) ~ "F"
  ))

# Factor
toFactor <- c('sex', 'study', 'scanner_id', 'Diagnosis')
enigmaDf[toFactor] <- lapply(enigmaDf[toFactor], factor)

# Add reason for scan column
clin22qFsDf$scan_reason_primary <- ''
clin22qIfsDf$scan_reason_primary <- ''
enigmaDf$scan_reason_primary <- ''

# Get the columns that are common to all dataframes
metadataCols <- c("patient_id", "age_at_scan_days", "age_in_years", 
                  "scan_id", "sex", "scanner_id", 
                  "Diagnosis", "FSVersion", 'scan_reason_primary')
globalPhenotypeCols <- c("BrainSegVol", "CerebralWhiteMatterVol",  "TotalGrayVol", 
                         "eTIV", "SurfaceHoles", "SubCortGrayVol",
                         "VentricleVolume", "CorticalSurfaceArea", 
                         "MeanCorticalThickness") #- missing in ifs?

clipFsDf <- getRegionalPhenotypes(clipFsDf, c(metadataCols, globalPhenotypeCols))
clipIfsDf <- getRegionalPhenotypes(clipIfsDf, c(metadataCols, globalPhenotypeCols))
clin22qFsDf <- getRegionalPhenotypes(clin22qFsDf, c(metadataCols, globalPhenotypeCols))
clin22qIfsDf <- getRegionalPhenotypes(clin22qIfsDf, c(metadataCols, globalPhenotypeCols))

# enigmaDf <- getRegionalPhenotypes(enigmaDf, c(metadataCols, globalPhenotypeCols))
enigmaDf <- enigmaDf[, c(metadataCols, globalPhenotypeCols)]

# Combine the data frames
myData <- rbind(clipFsDf, clipIfsDf, clin22qFsDf, clin22qIfsDf) #, enigmaDf)
myData$scan_reason_primary <- as.character(myData$scan_reason_primary)
myData[myData$scan_reason_primary == 'lesion', ]$scan_reason_primary <- 'non-brain lesion'


# Set the strings that were just added to be factors
myData$Diagnosis <- as.factor(myData$Diagnosis)
myData$FSVersion <- as.factor(myData$FSVersion)
myData$scanner_id <- as.factor(myData$scanner_id)

# Add a group column to both clinical and enigma datafranes
myData$Group <- paste(as.character(myData$Diagnosis), as.character(myData$FSVersion))
myData$Group <- as.factor(myData$Group)
enigmaDf$Group <- paste(as.character(enigmaDf$Diagnosis), as.character(enigmaDf$FSVersion))
enigmaDf$Group <- as.factor(enigmaDf$Group)

# Drop any rows missing sex and refactor the column
myData <- myData[myData$sex != '', ]
myData$sex <- as.factor(myData$sex)

# Sort the dataframe by patient_id and scanner_id
myData <- myData[ with(myData, order(myData$patient_id, myData$scan_id)), ]
myData <- myData[!duplicated(myData$patient_id), ]
myData$patient_id <- droplevels(myData$patient_id)

# maxAge <- max(enigmaDf$age_in_years, myData$age_in_years)
myData <- myData[ myData$age_at_scan_days > 55, ]
maxAge <- max(myData$age_in_years)

# Step 2: Generate distributions of age and sex for CLIP and 22q studies ----------------
png('/Users/youngjm/Data/22q11/figures/2022-05-09_demographics.png', width=165, height=135, units='mm', res=300)
# Create the age/sex plots
pAgesHC <- generateAgeSexDistributionPlot(myData, 'HC', myPalette, "Clinical", maxAge)
pAges22q <- generateAgeSexDistributionPlot(myData, '22q11DS', myPalette, "Clinical", maxAge)

# Combine and display the plots
grob <- patchworkGrob(pAgesHC + fig_pub_theme + pAges22q + fig_pub_theme + plot_layout(ncol = 1) + fig_pub_theme)
gridExtra::grid.arrange(grob)
dev.off()
# ggsave('/Users/youngjm/Data/22q11/figures/2022-05-09_demographics.png', width=7, height=7, units='in', device = 'png', dpi=2400)


# Make a plot for the top scan reasons 
myData$scan_reason_primary <- as.factor(myData$scan_reason_primary)
myData <- addPrimaryScanReasonCol(myData[myData$Diagnosis == 'HC', ])
reasonsTable <- summary(myData[myData$Diagnosis == 'HC', ]$top_scan_reason_factors)
pieLabels <- paste(names(reasonsTable), reasonsTable, sep=': ')
pie(reasonsTable, pieLabels, col=cbbPalette, mai=c(0,0,0,0),
    main="Top Scan Reasons for Clinical Controls (with sample sizes)")
ggsave('/Users/youngjm/Data/22q11/figures/2022-05-09_clip_reasons_for_scan.png', width=10, height=8.5, units='in', device = 'png', dpi=2400)

# Step 3: ComBat the CLIP and 22q data -----------------------------------------------------

phenotypes <- c("TotalGrayVol", "CerebralWhiteMatterVol", "SubCortGrayVol", "eTIV",
                "VentricleVolume", "CorticalSurfaceArea", "MeanCorticalThickness")

# Pull the scan ids
scanIds <- myData$scan_id
toCombat <- myData[phenotypes]
# Set the row names of the isolated phenotypes
toCombat$scan_id <- scanIds
# We want to remove differences based on the scanner id, so pull that info to use as the batch variable
# We want to protect for age, sex, and dx
covars <- data.frame(SITE = myData$scanner_id,
                     age_in_years = myData$age_in_years,
                     sex = myData$sex,
                     Diagnosis = myData$Diagnosis)

colnames(covars)
colnames(toCombat)

# Save the data and covars dataframes to csvs
write.csv(toCombat,"/Users/youngjm/Data/22q11/derivatives/mpr_fs_reconall_7.1.1_tables/2022-04_clip_22q_global_phenotypes_toCombat.csv", row.names = FALSE)
write.csv(covars,"/Users/youngjm/Data/22q11/derivatives/mpr_fs_reconall_7.1.1_tables/2022-04_clip_22q_covariates_toCombat.csv", row.names = FALSE)

combattedDf <- read.csv('/Users/youngjm/Data/22q11/derivatives/mpr_fs_reconall_7.1.1_tables/2022-04_clip_22q_global_phenotypes_combatted.csv')

# Add in other metadata needed for analysis
# combattedDf$scan_id <- scanIds
metaCols <- append(colnames(myData)[1:3], colnames(myData)[5:8])
metadataDf <- myData[, metaCols]
combattedDf <- cbind(metadataDf, combattedDf)
combattedDf$SurfaceHoles <- myData$SurfaceHoles
combattedDf$age_in_years <- myData$age_in_years
combattedDf$Group <- myData$Group

# Step 4: Generate GAMS for raw and combat CLIP and 22q ------------------------

# Build a table with all of the global phenotype p-vals parameters
tableColumns <- c("Total Gray\nVolume", "Cerebral White\nMatter Volume",
                  "Subcortical Gray\nVolume", "eTIV",
                  "Ventricle\nVolume",  "Total Cortical\nSurface Area",
                  "Mean Cortical\nThickness")
tableRows <- c("Surface Holes", "Diagnosis",
               "Sex = Male", "Age") #, "Age by Dx")

generateHorizontalTrajectoriesAndTables(myData, phenotypes, tableColumns, tableRows, 'Raw Phenotypes')
generateHorizontalTrajectoriesAndTables(combattedDf, phenotypes, tableColumns, tableRows, 'Combatted Phenotypes')

generateVerticalTrajectories(combattedDf, phenotypes, tableColumns, tableRows, 'Combatted Phenotypes')

generateTableScaledBetas(combattedDf, phenotypes, tableColumns, tableRows, 'Combatted Phenotypes')

# Step 5: Generate regional phenotype plots ------------------------------------

hcOut <- generateHemispherePlots(myData[myData$Diagnosis == 'HC', ], 'Clinical Healthy Control')
q22Out <- generateHemispherePlots(myData[myData$Diagnosis == '22q11DS', ], 'Clinical 22q11DS')

hcRegionalPeaks <- hcOut[[1]]
q22RegionalPeaks <- q22Out[[1]]

hcRegionalPlot <- hcOut[[2]]
q22RegionalPlot <- q22Out[[2]]

hcRegionalPeaks$peaksDiff <- abs(as.numeric(hcRegionalPeaks$em) - as.numeric(q22RegionalPeaks$em))
# Plot the dataframe using ggseg
diffPlot <- hcRegionalPeaks %>% 
  ggseg(mapping=aes(fill=as.numeric(peaksDiff)),
        hemisphere='left') +
  labs(title = "Diff. Age at Peak (years)", legend="Age") +
  theme(axis.title = element_blank()) +
  scale_fill_gradient(low = "blue", high = "red", na.value = NA) 

tmp <- list(hcRegionalPlot, q22RegionalPlot, diffPlot)
p <- wrap_plots(tmp, nrow=3, guides = "collect")
p

# Step 6: Generate distributions of age and sex for ENIGMA  --------------------

# Create the age/sex plots
pAgesEnigmaHC <- generateAgeSexDistributionPlot(enigmaDf, 'HC', myPalette, 'ENIGMA', maxAge)
pAgesEnigma22q <- generateAgeSexDistributionPlot(enigmaDf, '22q11DS', myPalette, 'ENIGMA', maxAge)

# Combine and display the plots
grob <- patchworkGrob(pAgesEnigmaHC + pAgesEnigma22q + plot_layout(guides="collect", ncol = 2))
gridExtra::grid.arrange(grob)

# Step 7: ComBat the ENIGMA data -----------------------------------------------

# Pull the scan ids
eScanIds <- enigmaDf$scan_id
toCombat <- enigmaDf[phenotypes]
# Set the row names of the isolated phenotypes
toCombat$scan_id <- eScanIds
# We want to remove differences based on the scanner id, so pull that info to use as the batch variable
# We want to protect for age, sex, and dx
covars <- data.frame(SITE = enigmaDf$scanner_id,
                     age_in_years = enigmaDf$age_in_years,
                     sex = enigmaDf$sex,
                     Diagnosis = enigmaDf$Diagnosis)

colnames(covars)
colnames(toCombat)

# Save the data and covars dataframes to csvs
# write.csv(toCombat,"/Users/youngjm/Data/22q11/tables/2022-04_enigma_global_phenotypes_toCombat.csv", row.names = FALSE)
# write.csv(covars,"/Users/youngjm/Data/22q11/tables/2022-04_enigma_covariates_toCombat.csv", row.names = FALSE)

combattedEnigmaDf <- read.csv('/Users/youngjm/Data/22q11/tables/2022-04_enigma_global_phenotypes_combatted.csv')

# Add in other metadata needed for analysis
combattedEnigmaDf$scan_id <- eScanIds
metaCols <- c(colnames(enigmaDf)[1:3], colnames(enigmaDf)[5:9])
metadataDf <- enigmaDf[, metaCols]
combattedEnigmaDf <- cbind(metadataDf, combattedEnigmaDf)
combattedEnigmaDf$SurfaceHoles <- enigmaDf$SurfaceHoles
combattedEnigmaDf$age_in_years <- enigmaDf$age_in_years
combattedEnigmaDf$Group <- enigmaDf$Group

# Step 8: Generate GAMS for raw and combat ENIGMA ------------------------------
generateHorizontalTrajectoriesAndTables(enigmaDf, phenotypes, tableColumns, tableRows, 'Raw ENIGMA Phenotypes')
generateHorizontalTrajectoriesAndTables(combattedEnigmaDf, phenotypes, tableColumns, tableRows, 'Combatted ENIGMA Phenotypes')

# Step 9: Combine ENIGMA, CLIP, and 22q ----------------------------------------
clipFsDf <- clipFsDf[, c(metadataCols, globalPhenotypeCols)]
clipIfsDf <- clipIfsDf[,  c(metadataCols, globalPhenotypeCols)]
clin22qFsDf <- clin22qFsDf[, c(metadataCols, globalPhenotypeCols)]
clin22qIfsDf <- clin22qIfsDf[, c(metadataCols, globalPhenotypeCols)]
enigmaDf <- enigmaDf[, c(metadataCols, globalPhenotypeCols)]

allData <- rbind(clipFsDf, clipIfsDf, clin22qFsDf, clin22qIfsDf, enigmaDf)
allData$Group <- paste(as.character(allData$Diagnosis), as.character(allData$FSVersion))
allData$Group <- as.factor(allData$Group)
allData$Diagnosis <- as.factor(allData$Diagnosis)
allData$FSVersion <- as.factor(allData$FSVersion)

# Step 10: Generate distributions of age and sex for all ------------------------

# Create the age/sex plots
pAgesAllHC <- generateAgeSexDistributionPlot(allData, 'HC', myPalette, "Clinical + ENIGMA", maxAge)
pAgesAll22q <- generateAgeSexDistributionPlot(allData, '22q11DS', myPalette, "Clinical + ENIGMA", maxAge)

# Combine and display the plots
grob <- patchworkGrob(pAgesHC + pAges22q + 
                      pAgesEnigmaHC + pAgesEnigma22q + 
                      pAgesAllHC + pAgesAll22q + plot_layout(guides="collect", ncol = 2))
gridExtra::grid.arrange(grob)

# Step 11: ComBat all ----------------------------------------------------------

# Pull the scan ids
aScanIds <- allData$scan_id
toCombat <- allData[phenotypes]
# Set the row names of the isolated phenotypes
toCombat$scan_id <- aScanIds
# We want to remove differences based on the scanner id, so pull that info to use as the batch variable
# We want to protect for age, sex, and dx
covars <- data.frame(SITE = allData$scanner_id,
                     age_in_years = allData$age_in_years,
                     sex = allData$sex,
                     Diagnosis = allData$Diagnosis)

colnames(covars)
colnames(toCombat)

# Save the data and covars dataframes to csvs
write.csv(toCombat,"/Users/youngjm/Data/22q11/tables/2022-04_clip_22q_enigma_global_phenotypes_toCombat.csv", row.names = FALSE)
write.csv(covars,"/Users/youngjm/Data/22q11/tables/2022-04_clip_22q_enigma_covariates_toCombat.csv", row.names = FALSE)

combattedAllDf <- read.csv('/Users/youngjm/Data/22q11/tables/2022-04_clip_22q_enigma_global_phenotypes_combatted.csv')

# Add in other metadata needed for analysis
combattedAllDf$scan_id <- aScanIds
metaCols <- append(colnames(allData)[1:3], colnames(allData)[5:8])
metadataDf <- allData[, metaCols]
combattedAllDf <- cbind(metadataDf, combattedAllDf)
combattedAllDf$SurfaceHoles <- allData$SurfaceHoles
combattedAllDf$age_in_years <- allData$age_in_years
combattedAllDf$Group <- allData$Group

# Step 12: Generate GAMs for raw and combat all --------------------------------
generateHorizontalTrajectoriesAndTables(allData, phenotypes, tableColumns, tableRows, 'All Raw Phenotypes')
generateHorizontalTrajectoriesAndTables(combattedAllDf, phenotypes, tableColumns, tableRows, 'All Combatted Phenotypes')


#-------------------------------------------------------------------------------
