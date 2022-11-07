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
# Loading and Prepping Data
#-------------------------------------------------------------------------------



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

tableRows <- c("Surface Holes", 
               "Age")
combatPlotsTablesM <- generatePlotsAndTablesForDataset(combattedDf, miniGlobalCols, tableColumns, tableRows, 'Combatted Phenotypes', "M")
combatPlotsTablesF <- generatePlotsAndTablesForDataset(combattedDf, miniGlobalCols, tableColumns, tableRows, 'Combatted Phenotypes', "F")
# combatHighQPlotsTables <- generatePlotsAndTablesForDataset(combattedHighQDf, miniGlobalCols, tableColumns, tableRows, 'Combatted Phenotypes')
# combatSuperHighQPlotsTables <- generatePlotsAndTablesForDataset(combattedSuperHighQDf, miniGlobalCols, tableColumns, tableRows, 'Combatted Phenotypes')

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

# Calculate the correlation between Lifespan and CLIP age at peak 
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





## Step ???: PCA? -------------

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
# predCentFn <- "/Users/youngjm/Data/lifespan_growth_charts/lifespan_centile_medians.csv"
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

modelFixedValues <- list(SurfaceHoles = mean(combattedDf$SurfaceHoles), 
                         sex='M', 
                         top_scan_reason_factors='headaches')

regionalTrajectoryPlots <- c()

for (i in 1:length(parsedRegionalPhenotypes)) {
  pheno <- parsedRegionalPhenotypes[i]
  print(pheno)
  combattedDf[[pheno]] <- (combattedDf[[paste0('lh_', pheno, '_grayVol')]] + 
    combattedDf[[paste0('rh_', pheno, '_grayVol')]])/2.0
  
  gam <- createMainGamm(combattedDf, pheno)
  preds <- predict_gam(gam, values = modelFixedValues)
  plotScatterScanReason <- generatePlotScatterWithCI(preds, combattedDf, 
                                                     combattedDf[ , pheno], pheno)
  regionalTrajectoryPlots[[pheno]] <- plotScatterScanReason
}


patchwork <- wrap_plots(regionalTrajectoryPlots, byrow=TRUE, guides="collect")
patchwork

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
