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
                   'sub-22q0295_ses-3884age02782_acq-MPROriginal3p0UnknownContrastHighResFromScanner45428_run-02_T1w',
                   'sub-22q0218_ses-4226age00675_acq-MPROriginal3p0UnknownContrastHighResFromScanner40160_run-01_T1w')
#------------------------------------------------------------------------------
# Step 02: combine the data frames into a single data frame
#------------------------------------------------------------------------------

# Add a Group column to each data frame to specify the group the data belongs to
clip.fs.data$Group <- 'Control FS'
clip.ifs.data$Group <- 'Control IFS'
q22.fs.data$Group <- '22qDS FS'
q22.ifs.data$Group <- '22qDS IFS'

# Add a Diagnosis column to each data frame to specify the group the data belongs to
clip.fs.data$Diagnosis <- 'Control'
clip.ifs.data$Diagnosis <- 'Control'
q22.fs.data$Diagnosis <- '22qDS'
q22.ifs.data$Diagnosis <- '22qDS'

# Add a Processing column to each data frame to specify the group the data belongs to
clip.fs.data$Processing <- 'FS'
clip.ifs.data$Processing <- 'IFS'
q22.fs.data$Processing <- 'FS'
q22.ifs.data$Processing <- 'IFS'

labels <- c('22qDS FS', '22qDS IFS', 'Control FS', 'Control IFS')

# Combine the data frames
myData <- rbind(clip.fs.data, clip.ifs.data, q22.fs.data, q22.ifs.data)

# Set the strings that were just added to be factors
myData$Group <- as.factor(myData$Group)
myData$Diagnosis <- as.factor(myData$Diagnosis)
myData$Processing <- as.factor(myData$Processing)

# Add a column with the log of age
myData$log_age_at_scan_days <- log(myData$age_at_scan_days)

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
    Diagnosis == '22qDS' ~ 1
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

# Generate GAMMs
gammTGV <- createGamm(myData, "TotalGray")
gammBrainSeg <- createGamm(myData, "BrainSeg")
gammCWM <- createGamm(myData, "CerebralWhiteMatter")
gammETIV <- createGamm(myData, "EstimatedTotalIntraCranialVol")
gammCSF <- createGamm(myData, "CSF")
gammSubCortGray <- createGamm(myData, "SubCortGrayVol")
gammCortSurf <- createGamm(myData, "SumCorticalSurfaceArea")
gammCortThick <- createGamm(myData, "AvgCorticalThickAvg")

# Summarize GAMMS
# summary(gammTGV$gam)
# summary(gammBrainSeg$gam)
# summary(gammCWM$gam)
# summary(gammETIV$gam)
# summary(gammCSF$gam)
# summary(gammSubCortGray$gam)
# summary(gammCortSurf$gam)
# summary(gammCortThick$gam)

# Predict on models
modelFixedValues <- list(SurfaceHoles = mean(myData$SurfaceHoles), sex='M', scanner_id='35008')
gammTGVPreds <- predict_gam(gammTGV$gam, values = modelFixedValues) #, scanner_id=names(sort(summary(myData$scanner_id), decreasing=T)[1])))# + resid(test)
gammBrainSegPreds <- predict_gam(gammBrainSeg$gam, values = modelFixedValues) #, scanner_id=names(sort(summary(myData$scanner_id), decreasing=T)[1])))# + resid(test)
gammCWMPreds <- predict_gam(gammCWM$gam, values = modelFixedValues) #, scanner_id=names(sort(summary(myData$scanner_id), decreasing=T)[1])))# + resid(test)
gammETIVPreds <- predict_gam(gammETIV$gam, values = modelFixedValues) #, scanner_id=names(sort(summary(myData$scanner_id), decreasing=T)[1])))# + resid(test)
gammCSFPreds <- predict_gam(gammCSF$gam, values = modelFixedValues)
gammSubCortGrayPreds <- predict_gam(gammSubCortGray$gam, values = modelFixedValues)
gammCortSurfPreds <- predict_gam(gammCortSurf$gam, values = modelFixedValues)
gammCortThickPreds <- predict_gam(gammCortThick$gam, values = modelFixedValues)

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
       fill="Sex & Diagnosis")

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
    # scale_x_continuous(trans='log10') +
    facet_wrap(~Diagnosis) +
    theme(axis.title = element_blank()) +
    labs(title = measureTitle)
  
  return(plot01)
}


generateOHBMPrettyPlotCI <- function(predictions, origData, measure, measureTitle){
  plot02 <- predictions %>%
    ggplot(aes(x=age_at_scan_days, y=fit, color=Diagnosis), log="x") +
    geom_smooth_ci() +
    scale_color_manual(values = dxColors) +
    scale_x_continuous(trans='log10') +
    theme(axis.title = element_blank()) +
    labs(title = "Growth Trajectory (log(Age))") #paste(measureTitle, ": Growth Trajectories"))
  
  return(plot02)
}

generateLinearTable <- function(aGam, nTests){
  subTab <- signif(summary(aGam)$p.table, 3)
  rows <- c('ordered(Diagnosis).L')
  cols <- c('t value', 'Pr(>|t|)')
  subTab <- subTab[,cols]
  subTab <- subset(subTab, rownames(subTab) %in% rows)
  # Add the p.adjust
  subTab <- cbind(subTab, p.adjust(subTab[,'Pr(>|t|)'], n = nTests))
  linT <- gridExtra::tableGrob(subTab, 
                               rows = 'Diagnosis',
                               cols = c('t value', 'p-value', 'Adj. p-value'))
  return(linT)
}

generateParametricTable <- function(aGam, nTests){
  subTab <- signif(summary(aGam)$s.table, 3)
  cols <- c("edf", "F","p-value")
  subTab <- subTab[,cols]
  rows <- c("s(log(age_at_scan_days)):ordered(Diagnosis)Control")
  subTab <- subset(subTab, rownames(subTab) %in% rows)
  # Add the p.adjust
  subTab <- cbind(subTab, p.adjust(subTab[,"p-value"], n = nTests))
  paraT <- gridExtra::tableGrob(subTab,
                                rows = 'Age by Diagnosis',
                                cols = c("edf", "F","p-value", 'Adj. p-value'))
}

# Plot pretty model graphs
nTests <- 6
p1 <- generateOHBMPrettyPlotScatter(gammTGVPreds, myData, myData$TotalGray, 'Total Gray Volume')
p2 <- generateOHBMPrettyPlotCI(gammTGVPreds, myData, myData$TotalGray, 'Total Gray Volume')
t1 <- generateParametricTable(gammTGV$gam, nTests)
t2 <- generateLinearTable(gammTGV$gam, nTests)

p3 <- generateOHBMPrettyPlotScatter(gammSubCortGrayPreds, myData, myData$SubCortGrayVol, 'SubCortical Gray Volume')
p4 <- generateOHBMPrettyPlotCI(gammSubCortGrayPreds, myData, myData$SubCortGrayVol, 'SubCortical Gray Volume')
t3 <- generateParametricTable(gammSubCortGray$gam, nTests)
t4 <- generateLinearTable(gammSubCortGray$gam, nTests)

p5 <- generateOHBMPrettyPlotScatter(gammCWMPreds, myData, myData$CerebralWhiteMatter, 'Cerebral White Matter Volume')
p6 <- generateOHBMPrettyPlotCI(gammCWMPreds, myData, myData$CerebralWhiteMatter, 'Cerebral White Matter Volume')
t5 <- generateParametricTable(gammCWM$gam, nTests)
t6 <- generateLinearTable(gammCWM$gam, nTests)

p7 <- generateOHBMPrettyPlotScatter(gammCSFPreds, myData, myData$CSF, 'CSF Volume')
p8 <- generateOHBMPrettyPlotCI(gammCSFPreds, myData, myData$CSF, 'CSF Volume')
t7 <- generateParametricTable(gammCSF$gam, nTests)
t8 <- generateLinearTable(gammCSF$gam, nTests)

p9 <- generateOHBMPrettyPlotScatter(gammCortSurfPreds, myData, myData$SumCorticalSurfaceArea, 'Total Cortical Surface Area')
p10 <- generateOHBMPrettyPlotCI(gammCortSurfPreds, myData, myData$SumCorticalSurfaceArea, 'Total Cortical Surface Area')
t9 <- generateParametricTable(gammCortSurf$gam, nTests)
t10 <- generateLinearTable(gammCortSurf$gam, nTests)

p11 <- generateOHBMPrettyPlotScatter(gammCortThickPreds, myData, myData$AvgCorticalThickAvg, 'Average Cortical Thickness (All Regions)')
p12 <- generateOHBMPrettyPlotCI(gammCortThickPreds, myData, myData$AvgCorticalThickAvg, 'Average Cortical Thickness (All Regions)')
t11 <- generateParametricTable(gammCortThick$gam, nTests)
t12 <- generateLinearTable(gammCortThick$gam, nTests)


pgrobble <- patchworkGrob(p1 + p2 +
                           p3 + p4 + 
                           p5 + p6 + 
                           p7 + p8 + 
                           p9 + p10 + 
                           p11 + p12 + 
                           plot_layout(guides="collect", widths = c(2, 1), ncol = 2))

tg1 <- gtable_combine(t1,t2, along=2)
tg2 <- gtable_combine(t3,t4, along=2)
tg3 <- gtable_combine(t5,t6, along=2)
tg4 <- gtable_combine(t7,t8, along=2)
tg5 <- gtable_combine(t9,t10, along=2)
tg6 <- gtable_combine(t11,t12, along=2)

# title <- textGrob("Title",gp=gpar(fontsize=14))
# table <- gtable_add_rows(tg1, 
#                          heights = grobHeight(title))
# table <- gtable_add_grob(table, list(title, footnote),
#                          t=c(1, nrow(table)), l=c(1,2), 
#                          r=ncol(table))

lay <- rbind(c(1,2),
             c(1,3),
             c(1,4),
             c(1,5),
             c(1,6),
             c(1,7))

gridExtra::grid.arrange(pgrobble, t1, t3, t5, t7, t9, t11,
                        top = textGrob("Age vs. Phenotype", gp=gpar(fontsize=20)),
                        left = textGrob("Phenotype", gp=gpar(fontsize=18), rot=90), 
                        bottom = textGrob("Age at Scan (Days)", gp=gpar(fontsize=18)),
                        layout_matrix=lay, widths=c(3,1))

