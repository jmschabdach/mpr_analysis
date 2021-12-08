library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgcv)
library(tidymv)
library(patchwork) # graph organization within a figure
library(gtsummary)
library(grid)
# library(huxtable)
# library(magrittr)
# library(MatchIt)

#------------------------------------------------------------------------------
# Step 01: Read the data
#
# - Specify the file names
# - Read the csv files
# - Print the first 10 lines of each file
#------------------------------------------------------------------------------

# Specify the file names
clip.fs.fn <- '/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_7.1.1_structural_stats.csv'
clip.ifs.fn <- '/Users/youngjm/Data/clip/images/derivatives/mpr_ifs_reconall_7.1.1_structural_stats.csv'
q22.fs.fn <- '/Users/youngjm/Data/22q11/derivatives/mpr_fs_reconall_7.1.1_structural_stats.csv'
q22.ifs.fn <- '/Users/youngjm/Data/22q11/derivatives/mpr_ifs_reconall_7.1.1_structural_stats.csv'

# Read the data from the files
clip.fs.data <- read.csv(clip.fs.fn, stringsAsFactors = TRUE)
clip.ifs.data <- read.csv(clip.ifs.fn, stringsAsFactors = TRUE)
q22.fs.data <- read.csv(q22.fs.fn, stringsAsFactors = TRUE)
q22.ifs.data <- read.csv(q22.ifs.fn, stringsAsFactors = TRUE)

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
  
  return(new.data.df)
}

clip.fs.data <- prepInfo(clip.fs.data)
clip.ifs.data <- prepInfo(clip.ifs.data)
q22.fs.data <- prepInfo(q22.fs.data)
q22.ifs.data <- prepInfo(q22.ifs.data)

# Drop elements with scanIds

scanIdsToDrop <- c('sub-22q0091_ses-8236age11798_acq-MPROriginal1p5UnknownContrastHighResFromScanner7315_run-01_T1w',
                   'sub-22q0172_ses-5478age03490_acq-MPRDerived3p0PostcontrastHighResFromScanner67016_run-05_T1w',
                   'sub-22q0250_ses-2074age03002_acq-MPRDerived3p0PostcontrastHighResFromScanner40156_run-04_T1w',
                   'sub-22q0250_ses-2074age03002_acq-MPRDerived3p0PostcontrastHighResFromScanner40156_run-05_T1w',
                   'sub-22q0250_ses-2074age03002_acq-MPRDerived3p0PrecontrastHighResFromScanner40156_run-02_T1w',
                   'sub-22q0265_ses-2436age02562_acq-MPRDerived1p5UnknownContrastHighResFromScanner169573_run-02_T1w',
                   'sub-22q0298_ses-1529age00771_acq-MPRDerived1p5UnknownContrastHighResFromScanner169573_run-02_T1w',
                   'sub-22q0300_ses-2415age01053_acq-MPRDerived3p0UnknownContrastHighResFromScanner0000000CHOPPETMR_run-02_T1w',
                   'sub-22q0300_ses-2415age01053_acq-MPRDerived3p0UnknownContrastHighResFromScanner0000000CHOPPETMR_run-03_T1w',
                   'sub-22q0305_ses-1026age00874_acq-MPRDerived1p5UnknownContrastHighResFromScanner169573_run-02_T1w',
                   'sub-22q0313_ses-3555age02342_acq-MPRDerived3p0UnknownContrastHighResFromScanner0000000CHOPPETMR_run-02_T1w',
                   'sub-22q0313_ses-3555age02342_acq-MPRDerived3p0UnknownContrastHighResFromScanner0000000CHOPPETMR_run-03_T1w',
                   'sub-22q0126_ses-7492age06097_acq-MPROriginal3p0UnknownContrastHighResFromScanner35069_run-01_T1w', #only one on this scanner
                   'sub-22q0295_ses-3884age02782_acq-MPROriginal3p0UnknownContrastHighResFromScanner45428_run-01_T1w',
                   'sub-22q0218_ses-4226age00675_acq-MPROriginal3p0UnknownContrastHighResFromScanner40160_run-01_T1w')
#------------------------------------------------------------------------------
# Step 02: combine the data frames into a single data frame
#------------------------------------------------------------------------------

# Add a Group column to each data frame to specify the group the data belongs to
clip.fs.data$Group <- 'CLIP FS'
clip.ifs.data$Group <- 'CLIP IFS'
q22.fs.data$Group <- '22q FS'
q22.ifs.data$Group <- '22q IFS'

# Add a Diagnosis column to each data frame to specify the group the data belongs to
clip.fs.data$Diagnosis <- 'CLIP'
clip.ifs.data$Diagnosis <- 'CLIP'
q22.fs.data$Diagnosis <- '22q'
q22.ifs.data$Diagnosis <- '22q'

# Add a Processing column to each data frame to specify the group the data belongs to
clip.fs.data$Processing <- 'FS'
clip.ifs.data$Processing <- 'IFS'
q22.fs.data$Processing <- 'FS'
q22.ifs.data$Processing <- 'IFS'

# vivid blue for the CLIP FS
# sky blue for the CLIP IFS
# orange for the 22q FS
# peach for the 22q IFS
my.colors <- c('#f97306', '#ffb07c', '#152eff', '#75bbfd')

labels <- c('22q FS', '22q IFS', 'CLIP FS', 'CLIP IFS')

# Combine the data frames
myData <- rbind(clip.fs.data, clip.ifs.data, q22.fs.data, q22.ifs.data)
# clipData <- rbind(clip.fs.data, clip.ifs.data)
# q22Data <- rbind(q22.fs.data, q22.ifs.data)

# Set the strings that were just added to be factors
myData$Group <- as.factor(myData$Group)
myData$Diagnosis <- as.factor(myData$Diagnosis)
myData$Processing <- as.factor(myData$Processing)

#------------------------------------------------------------------------------
# Step 02b: clean the data (only first scan of 3T, only one scan per subject, 
#           Diagnosis as 0/1)
#------------------------------------------------------------------------------

# Only the 3.0T scans
myData <- myData[myData$MagneticFieldStrength == '3.0', ]

# drop anything in the scanIdsToDrop
myData$scan_id <- as.character(myData$scan_id)
myData <- myData[!myData$scan_id %in% scanIdsToDrop, ]

# Make a Has22q column
myData <- myData %>%
  mutate(Has22q = case_when(
    Diagnosis == 'CLIP' ~ 0,
    Diagnosis == '22q' ~ 1
  ))

# Cross sectionalize the data
# First, order by subject and age
myData <- myData[with(myData, order(patient_id, age_at_scan_days, scan_id)), ]
# Now get first scan per subject
myData <- myData[!duplicated(myData$patient_id, myData$age_at_scan_days), ]

# After checking, found only 1 subject with unknown sex, so removing it
myData <- myData[myData$sex != 'U', ]
myData <- droplevels(myData)

#------------------------------------------------------------------------------
# Step 03: Replicate the Python plots
#
# - Plot age vs. eTIV as a scatter plot
# - Plot age vs. SurfaceHoles as a scatter plot
# - Plot histogram of SurfaceHoles as a QA measure
#------------------------------------------------------------------------------

# Plot age vs. eTIV as a scatter plot
ggplot(data=myData, aes(x=age_at_scan_days, y=EstimatedTotalIntraCranialVol, color = Group)) + # Set up the x and y data
  geom_point() + 
  scale_color_manual(values = my.colors) +
  labs(title = "Age vs. eTIV", 
       x = 'Age at Scan (Days)',
       y = 'Estimated Total Intracranial Volume')

# Not sure how to add a legend, but ...

# Plot age vs. SurfaceHoles as a scatter plot
ggplot(data=myData, aes(x=age_at_scan_days, y=SurfaceHoles, color = Group)) + # Set up the x and y data
  geom_point() + # add the scatter plot
  scale_color_manual(values = my.colors) + 
  labs(title = "Age vs. # Surface Holes", 
       x = 'Age at Scan (Days)',
       y = '# Surface Holes')


# Plot histogram of SurfaceHoles
ggplot(data=myData, aes(x=SurfaceHoles, fill=Group)) +
  geom_histogram(binwidth = 5) +
  geom_density(alpha=0.2) +
  scale_fill_manual(values = my.colors) +
  labs(title = "Histogram of # Surface Holes",
       x = '# Surface Holes',
       y = '# Occurrences')

#------------------------------------------------------------------------------
# Step 04: Basic linear model
#------------------------------------------------------------------------------

# Age, scanner, FS/IFS, sex, surface holes, diagnosistic group, age2, age3

## Create a basic linear model with a prespecified set of covariates
# @param df A data frame containing at least the columns 'age_at_scan_days' and 'SurfaceHoles'
# @return basicModel A List object produced by the lm function
createBasicLinearModel <- function(df, measure) {
  # Make a basic linear model accounting for age and surface holes for the given data frame
  basicModel <- lm(measure~age_at_scan_days+
                     (age_at_scan_days**2)+(age_at_scan_days**3)+
                     SurfaceHoles,
                   data = df)
  # Return the model
  return(basicModel)
}

clipGMModel <- createBasicLinearModel(myData, myData$TotalGray)

# Make a new dataframe predDf consisting of clipData with the three outputs of predict as new columns
# Predict with interval='confidence' produces 3 outputs: fit, lwr, and upr
predDf <- cbind(myData, predict(clipGMModel, interval='confidence'))

# Plot age vs. gray matter as a scatter plot
ggplot(data=predDf, aes(x=age_at_scan_days, y=TotalGray, color=Group)) + # Set up the default data frame and x and y data a
  geom_point() + # add the scatter plot, inherits the aes() from the ggplot
  geom_line(aes(x=age_at_scan_days, y=fit)) + # add the regression line, specify the x and the predicted y using aes()
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.3) + # add a partially shaded ribbon displaying the confidence interval
  scale_color_manual(values = my.colors) +
  facet_wrap(~Has22q+sex) +
  labs(title = "Age vs. Total Gray Volume",
       x = 'Age at Scan (Days)',
       y = 'Total Gray Volume')

#------------------------------------------------------------------------------
# Step 05: General Additive Model (GAM)
#------------------------------------------------------------------------------

createMultivariateGam <- function(df, measure) {
  # Make a basic linear model accounting for age and surface holes for the given data frame
  # s() smooths the covariate data
  # s() can take argument by, which groups by another covariate
  basicModel <- gam(formula=measure~s(log(age_at_scan_days), by = ordered(Diagnosis)) + # Using age gives squiggles, but log(age) doesn't
                      s(log(age_at_scan_days)) +
                      s(SurfaceHoles) + # definitely has an impact
                      ordered(Diagnosis) +
                      # scanner_id +
                      sex, # categorical = linear, don't smooth
                      # Diagnosis, # categorical = linear, don't smooth
                   data = df, method = "REML")
  # Return the model
  return(basicModel)
}

gam1 <- createMultivariateGam(myData, myData$TotalGray)
summary(gam1)

gamPredictions1 <- predict_gam(gam1, values = list(SurfaceHoles = mean(myData$SurfaceHoles), sex='M'))# + resid(test)

gamPredictions1 %>%
ggplot(aes(x=age_at_scan_days, y=fit)) +
  geom_point(data=myData, aes(x=age_at_scan_days, y=TotalGray, color=Group)) +
  geom_smooth_ci() +
  scale_color_manual(values = my.colors) +
  # facet_grid(sex ~ Diagnosis) +
  facet_wrap(~Diagnosis) +
  labs(title = "Age vs. Total Gray Volume",
       x = 'Age at Scan (Days)',
       y = 'Total Gray Volume')


#------------------------------------------------------------------------------
# Step 06: General Additive Mixed Model (GAMM)
#
# Slightly more complexity than a GAM. Will allow us to model the impact of
# individual scanners as random effects.
#------------------------------------------------------------------------------

createGamm <- function(df, measure) {
  formula <- as.formula(paste(measure, "s(log(age_at_scan_days), by = ordered(Diagnosis), fx=T) +
                      s(log(age_at_scan_days), fx=T) +
                      s(SurfaceHoles) +
                      ordered(Diagnosis) +
                      sex", sep="~"))
  
  # Make a basic linear model accounting for age and surface holes for the given data frame
  # s() smooths the covariate data
  # s() can take argument by, which groups by another covariate
  mixedModel <- gamm(formula,
                     random = list(scanner_id=~1),
                     data = df)
  # Return the model
  return(mixedModel)
}


generatePrettyPlot <- function(predictions, origData, measure, measureTitle){
  plot01 <- predictions %>%
    ggplot(aes(x=age_at_scan_days, y=fit)) +
    geom_point(data=origData, aes(x=age_at_scan_days, y=measure, color=Group)) +
    geom_smooth_ci() +
    scale_color_manual(values = my.colors) +
    facet_wrap(~Diagnosis) +
    theme(axis.title = element_blank()) +
    labs(title = "Models, CI, and Data")
  
  plot02 <- predictions %>%
    ggplot(aes(x=age_at_scan_days, y=fit, color=Diagnosis)) +
    geom_smooth_ci() +
    # scale_fill_manual(values = c('#f97306', '#152eff')) +
    theme(axis.title = element_blank()) +
    labs(title = "Models and CI")
  
  grob <- patchworkGrob(plot01 / plot02 + plot_layout(guides = 'collect', height = c(1, 2)))
  gridExtra::grid.arrange(grob, top = paste("Age vs.", measureTitle, sep=' '), left = measureTitle, bottom = "Age at Scan (Days)")
  
}

# Generate GAMMs
gammTGV <- createGamm(myData, "TotalGray")
gammBrainSeg <- createGamm(myData, "BrainSeg")
gammCWM <- createGamm(myData, "CerebralWhiteMatter")
gammETIV <- createGamm(myData, "EstimatedTotalIntraCranialVol")
# add CSF, subcortical gray, total surface area and mean cortical thickness
gammCSF <- createGamm(myData, "CSF")
gammSubCortGray <- createGamm(myData, "SubCortGrayVol")
gammCortSurf <- createGamm(myData, "SumCorticalSurfaceArea")
gammCortThick <- createGamm(myData, "AvgCorticalThickAvg")

# Summarize GAMMS
summary(gammTGV$gam)
summary(gammBrainSeg$gam)
summary(gammCWM$gam)
summary(gammETIV$gam)
summary(gammCSF$gam)
summary(gammSubCortGray$gam)
summary(gammCortSurf$gam)
summary(gammCortThick$gam)

# Predict on models
modelFixedValues <- list(SurfaceHoles = mean(myData$SurfaceHoles), sex='M')
gammTGVPreds <- predict_gam(gammTGV$gam, values = modelFixedValues) #, scanner_id=names(sort(summary(myData$scanner_id), decreasing=T)[1])))# + resid(test)
gammBrainSegPreds <- predict_gam(gammBrainSeg$gam, values = modelFixedValues) #, scanner_id=names(sort(summary(myData$scanner_id), decreasing=T)[1])))# + resid(test)
gammCWMPreds <- predict_gam(gammCWM$gam, values = modelFixedValues) #, scanner_id=names(sort(summary(myData$scanner_id), decreasing=T)[1])))# + resid(test)
gammETIVPreds <- predict_gam(gammETIV$gam, values = modelFixedValues) #, scanner_id=names(sort(summary(myData$scanner_id), decreasing=T)[1])))# + resid(test)
gammCSFPreds <- predict_gam(gammCSF$gam, values = modelFixedValues)
gammSubCortGrayPreds <- predict_gam(gammSubCortGray$gam, values = modelFixedValues)
gammCortSurfPreds <- predict_gam(gammCortSurf$gam, values = modelFixedValues)
gammCortThickPreds <- predict_gam(gammCortThick$gam, values = modelFixedValues)

# Plot pretty model graphs
generatePrettyPlot(gammTGVPreds, myData, myData$TotalGray, 'Total Gray Volume')
generatePrettyPlot(gammBrainSegPreds, myData, myData$BrainSeg, 'Brain Segmentation')
generatePrettyPlot(gammCWMPreds, myData, myData$CerebralWhiteMatter, 'Cerebral White Matter')
generatePrettyPlot(gammETIVPreds, myData, myData$EstimatedTotalIntraCranialVol, 'Estimated Total Intracranial Volume')
generatePrettyPlot(gammCSFPreds, myData, myData$CSF, 'CSF Volume')
generatePrettyPlot(gammSubCortGrayPreds, myData, myData$SubCortGrayVol, 'SubCortical Gray Volume')
generatePrettyPlot(gammCortSurfPreds, myData, myData$SumCorticalSurfaceArea, 'Total Cortical Surface Area')
generatePrettyPlot(gammCortThickPreds, myData, myData$AvgCorticalThickAvg, 'Average Cortical Thickness (All Regions)')


#-------------------------------------------------------------------------------
# Testing: identifying outliers
#-------------------------------------------------------------------------------

normalizeValue <- function(origValues){
  # denom
  denom = max(origValues) - min(origValues)
  # normalize the values
  for (i in 1:length(origValues)){
    origValues[[i]] <- (origValues[[i]] - min(origValues))/denom
  }
  
  return(origValues)
}

myData$normTGV <- normalizeValue(myData$TotalGray)
myData$normBrainSeg <- normalizeValue(myData$BrainSeg)
myData$normCWM <- normalizeValue(myData$CerebralWhiteMatter)
myData$normETIV <- normalizeValue(myData$EstimatedTotalIntraCranialVol)
myData$normCSF <- normalizeValue(myData$CSF)
myData$normSubCortGray <- normalizeValue(myData$SubCortGrayVol)
myData$normCortSurf <- normalizeValue(myData$SumCorticalSurfaceArea)
myData$normCortThick <- normalizeValue(myData$AvgCorticalThickAvg)

summedValues <- myData$normTGV + myData$normBrainSeg + myData$normCWM +
                myData$normETIV + myData$normCSF + myData$normSubCortGray +
                myData$normCortSurf + myData$normCortThick
summedValues <- summedValues/8
myData$trueNormal <- summedValues

# plot(myData$age_at_scan_days, myData$trueNormal)

myData[which(myData$age_at_scan_days > 750 & myData$trueNormal < 0.37),]
myData[myData$trueNormal > 0.9, ]

#-------------------------------------------------------------------------------
# Figures for OHBM Abstract
#-------------------------------------------------------------------------------


# Figure 1: 

myData$sexDxGroup <- as.factor(paste(as.character(myData$sex), as.character(myData$Diagnosis)))
sexColors <- c("#ffb2d0", "#ff69af", "#95d0fc", "#448ee4" )
dxColors <- c("#e66100", "#5d3a9b")

ages <- ggplot(data=myData, aes(x=age_at_scan_days, fill=sexDxGroup)) +
  geom_histogram(position="stack", binwidth=91) +
  scale_fill_manual(values = sexColors) +
  labs(title = "Distribution of Patient Age at Scan",
       x = 'Patient Age (Days)',
       y = '# Patients',
       fill="Sex + Dx")

scanners <- ggplot(data=myData, aes(x=scanner_id, fill=Diagnosis)) +
  geom_bar(position="stack") +
  scale_fill_manual(values = dxColors) +
  labs(title = "Distribution of Patients Between Scanners",
       x = "Scanner ID",
       y = "# Scans per Scanner")

grob <- patchworkGrob(ages + scanners + plot_layout(guides = 'collect'))
gridExtra::grid.arrange(grob)


generateOHBMPrettyPlotScatter <- function(predictions, origData, measure, measureTitle){
  plot01 <- predictions %>%
    ggplot(aes(x=age_at_scan_days, y=fit)) +
    geom_point(data=origData, aes(x=age_at_scan_days, y=measure, color=Diagnosis)) +
    geom_smooth_ci() +
    scale_color_manual(values = dxColors) +
    facet_wrap(~Diagnosis) +
    theme(axis.title = element_blank()) +
    labs(title = measureTitle)
  
  return(plot01)
}

generateOHBMPrettyPlotCI <- function(predictions, origData, measure, measureTitle){
  plot02 <- predictions %>%
    ggplot(aes(x=age_at_scan_days, y=fit, color=Diagnosis)) +
    geom_smooth_ci() +
    scale_color_manual(values = dxColors) +
    theme(axis.title = element_blank()) +
    labs(title = "Growth Trajectory") #paste(measureTitle, ": Growth Trajectories"))
  
  return(plot02)
}

# Plot pretty model graphs
p1 <- generateOHBMPrettyPlotScatter(gammTGVPreds, myData, myData$TotalGray, 'Total Gray Volume')
p2 <- generateOHBMPrettyPlotCI(gammTGVPreds, myData, myData$TotalGray, 'Total Gray Volume')

p3 <- generateOHBMPrettyPlotScatter(gammSubCortGrayPreds, myData, myData$SubCortGrayVol, 'SubCortical Gray Volume')
p4 <- generateOHBMPrettyPlotCI(gammSubCortGrayPreds, myData, myData$SubCortGrayVol, 'SubCortical Gray Volume')

p5 <- generateOHBMPrettyPlotScatter(gammCWMPreds, myData, myData$CerebralWhiteMatter, 'Cerebral White Matter Volume')
p6 <- generateOHBMPrettyPlotCI(gammCWMPreds, myData, myData$CerebralWhiteMatter, 'Cerebral White Matter Volume')

p7 <- generateOHBMPrettyPlotScatter(gammCSFPreds, myData, myData$CSF, 'CSF Volume')
p8 <- generateOHBMPrettyPlotCI(gammCSFPreds, myData, myData$CSF, 'CSF Volume')

p9 <- generateOHBMPrettyPlotScatter(gammCortSurfPreds, myData, myData$SumCorticalSurfaceArea, 'Total Cortical Surface Area')
p10 <- generateOHBMPrettyPlotCI(gammCortSurfPreds, myData, myData$SumCorticalSurfaceArea, 'Total Cortical Surface Area')

p11 <- generateOHBMPrettyPlotScatter(gammCortThickPreds, myData, myData$AvgCorticalThickAvg, 'Average Cortical Thickness (All Regions)')
p12 <- generateOHBMPrettyPlotCI(gammCortThickPreds, myData, myData$AvgCorticalThickAvg, 'Average Cortical Thickness (All Regions)')


grobble <- patchworkGrob(p1 + p2 + p3 + p4 + 
                         p5 + p6 + p7 + p8 + 
                         p9 + p10 + p11 + p12 + plot_layout(guides="collect", widths = c(2,1), ncol = 2))
gridExtra::grid.arrange(grobble, top = textGrob("Age vs. Phenotype", gp=gpar(fontsize=20)),
                                 left = textGrob("Phenotype", gp=gpar(fontsize=18), rot=90), 
                                 bottom = textGrob("Age at Scan (Days)", gp=gpar(fontsize=18)))

