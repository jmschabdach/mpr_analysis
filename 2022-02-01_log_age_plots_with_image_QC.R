library(ggplot2)
library(ggpubr)
library(dplyr)
library(mgcv)
library(tidymv)
library(patchwork) # graph organization within a figure
library(gtsummary)
library(grid)
library(harrypotter)
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
qc.fn <- '/Users/youngjm/Data/2021-12_mpr_analysis_primary_df_clip_22q.csv'
# This file was generated in the script 2021-12_07...r and is a combination of 
#  the imaging phenotypes of all images from the 22q and CLIP data

# Read the data from the files
qc.data <- read.csv(qc.fn, stringsAsFactors = TRUE)

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

qc.data <- prepInfo(qc.data)

# Add a column with the log of age
qc.data$log_age_at_scan_days <- log(qc.data$age_at_scan_days)

#------------------------------------------------------------------------------
# Step 02: clean the data (only first scan of 3T, only one scan per subject, 
#           Diagnosis as 0/1)
#------------------------------------------------------------------------------

# Only the 3.0T scans
qc.data <- qc.data[qc.data$MagneticFieldStrength == '3.0', ]

# Remove postcontrast scans
qc.data <- qc.data[qc.data$rawdata_image_grade > -1, ]

# Cross sectionalize the data
# First, order by subject and age and image grade
qc.data <- qc.data[with(qc.data, order(patient_id, age_at_scan_days, scan_id, rawdata_image_grade)), ]
# Now get first scan per subject
qc.data <- qc.data[!duplicated(qc.data$patient_id, qc.data$age_at_scan_days), ]

# Filter out unusable data grade 0
myData <- qc.data
myData <- qc.data[qc.data$rawdata_image_grade > 0, ]
# myData <- myData[myData$age_at_scan_days > 10, ]

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
                      rawdata_image_grade +
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
  scale_x_discrete(labels=c("20593" = "Scanner\n01", "20618" = "Scanner\n02",
                            "35008" = "Scanner\n03", "35014" = "Scanner\n04",
                            "40156" = "Scanner\n05", "40180" = "Scanner\n06",
                            "45195" = "Scanner\n07", "45886" = "Scanner\n08",
                            "46005" = "Scanner\n09", "67016" = "Scanner\n10")) +
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
  rows <- c('ordered(Diagnosis).L', 'rawdata_image_grade')
  cols <- c('t value', 'Pr(>|t|)')
  subTab <- subTab[,cols]
  subTab <- subset(subTab, rownames(subTab) %in% rows)
  # Add the p.adjust
  subTab <- cbind(subTab, p.adjust(subTab[,'Pr(>|t|)'], n = nTests))
  linT <- gridExtra::tableGrob(subTab, 
                               rows = c('Diagnosis', 'QC Rating'),
                               cols = c('t value', 'p-value', 'Adj. p-value'))
  return(linT)
}

generateParametricTable <- function(aGam, nTests){
  subTab <- signif(summary(aGam)$s.table, 3)
  cols <- c("edf", "F","p-value")
  subTab <- subTab[,cols]
  # rows <- c("s(log(age_at_scan_days)):ordered(Diagnosis)Control")
  # subTab <- subset(subTab, rownames(subTab) %in% rows)
  # Add the p.adjust
  subTab <- cbind(subTab, p.adjust(subTab[,"p-value"], n = nTests))
  paraT <- gridExtra::tableGrob(subTab,
                                # rows = 'Age by Diagnosis',
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

gridExtra::grid.arrange(pgrobble, t2, t4, t6, t8, t10, t12,
                        top = textGrob("Age vs. Phenotype", gp=gpar(fontsize=20)),
                        left = textGrob("Phenotype", gp=gpar(fontsize=18), rot=90), 
                        bottom = textGrob("Age at Scan (Days)", gp=gpar(fontsize=18)),
                        layout_matrix=lay, widths=c(3,1))

############################################3
#
############################################

qc.data$rawdata_image_grade <- as.factor(qc.data$rawdata_image_grade)

v <- qc.data %>%
  ggplot(aes(x=rawdata_image_grade, y=SurfaceHoles)) +
  geom_violin(aes(fill=rawdata_image_grade)) +
  geom_boxplot(width=0.1)+ 
  labs(title = "Distribution of Surface Holes Count at Scan Across Image QC Ratings",
       x = 'Image QC Rating Group',
       y = 'Number of Surface Holes (count)',
       fill="Image QC Rating Group")
v

ks.test(qc.data[qc.data$rawdata_image_grade == 0, ]$SurfaceHoles, qc.data[qc.data$rawdata_image_grade == 1, ]$SurfaceHoles)
ks.test(qc.data[qc.data$rawdata_image_grade == 0, ]$SurfaceHoles, qc.data[qc.data$rawdata_image_grade == 2, ]$SurfaceHoles)
ks.test(qc.data[qc.data$rawdata_image_grade == 2, ]$SurfaceHoles, qc.data[qc.data$rawdata_image_grade == 1, ]$SurfaceHoles)

# cor.test(as.numeric(qc.data$rawdata_image_grade), qc.data$SurfaceHoles)
mylogit <- glm(rawdata_image_grade ~ SurfaceHoles, data=qc.data, family='binomial')
summary(mylogit)




v <- qc.data %>%
  ggplot(aes(x=rawdata_image_grade, y=age_at_scan_days)) +
  geom_violin(aes(fill=rawdata_image_grade)) +
  geom_boxplot(width=0.1) + 
  labs(title = "Distribution of Patient Age at Scan Across Image QC Ratings",
       x = 'Image QC Rating Group',
       y = 'Age at Scan (days)',
       fill="Image QC Rating Group")
v

ks.test(qc.data[qc.data$rawdata_image_grade == 0, ]$age_at_scan_days, qc.data[qc.data$rawdata_image_grade == 1, ]$age_at_scan_days)
ks.test(qc.data[qc.data$rawdata_image_grade == 0, ]$age_at_scan_days, qc.data[qc.data$rawdata_image_grade == 2, ]$age_at_scan_days)
ks.test(qc.data[qc.data$rawdata_image_grade == 2, ]$age_at_scan_days, qc.data[qc.data$rawdata_image_grade == 1, ]$age_at_scan_days)

mylogit <- glm(rawdata_image_grade ~ age_at_scan_days, data=qc.data, family='binomial')
summary(mylogit)


############################################################
# Looking at IFS QC
############################################################

ifs.data <- qc.data[qc.data$Processing == 'IFS', ]

v <- ifs.data %>%
  ggplot(aes(x=rawdata_image_grade, y=age_at_scan_days)) +
  geom_violin(aes(fill=rawdata_image_grade)) +
  geom_boxplot(width=0.1)+ 
  labs(title = "Distribution of Age at Scan Across Image QC Ratings (Age < 3 Years)",
       x = 'Image QC Rating Group',
       y = 'Age at Scan (days)',
       fill="Image QC Rating Group")
v

mylogit <- glm(rawdata_image_grade ~ age_at_scan_days, data=ifs.data, family='binomial')
summary(mylogit)


##########################################################
# Plotting distribution of scan qc over age and scanner id
##########################################################

# sexColors <- c("#ffb2d0", "#ff69af", "#95d0fc", "#448ee4" )
# dxColors <- c("#e66100", "#5d3a9b")
qcColors <- c("#ffff66", "#c0c0c0", "#3333ff")

ages <- ggplot(data=qc.data, aes(x=age_at_scan_days, fill=as.factor(rawdata_image_grade))) +
  geom_histogram(position="stack", binwidth=91) +
  scale_fill_hp_d(option="Ravenclaw", name='Image Grade', direction=-1) +
  labs(title = "Distribution of Patient Age at Scan",
       x = 'Patient Age (Days)',
       y = '# Patients',
       fill="Image Grade")

scanners <- ggplot(data=qc.data, aes(x=scanner_id, fill=as.factor(rawdata_image_grade))) +
  geom_bar(position="stack") +
  scale_fill_hp_d(option="Ravenclaw", name='Image Grade', direction=-1) +
  scale_x_discrete(labels=c("20593" = "Scanner\n01", "20618" = "Scanner\n02",
                            "35008" = "Scanner\n03", "35014" = "Scanner\n04",
                            "40156" = "Scanner\n05", "40180" = "Scanner\n06",
                            "45195" = "Scanner\n07", "45886" = "Scanner\n08",
                            "46005" = "Scanner\n09", "67016" = "Scanner\n10")) +
  labs(title = "Distribution of Patients Between Scanners",
       x = "Scanner ID",
       y = "# Scans per Scanner",
       fill="Image Grade")

grob <- patchworkGrob(ages + scanners + plot_layout(guides = 'collect'))
gridExtra::grid.arrange(grob)
