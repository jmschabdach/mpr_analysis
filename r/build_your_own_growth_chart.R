#-------------------------------------------------------------------------------
# This script was assembled to walk a user through building a GAMLSS model
# using their own neuroimaging phenotype and demographic data. If you customize
# it and use it in your research, we politely ask that you cite the following 
# papers (at least the first one if your citations are limited):
#
# [SLIP paper](https://www.medrxiv.org/content/10.1101/2023.01.13.23284533v1)
# 
# [Lifespan Nature paper](https://www.nature.com/articles/s41586-022-04554-y)
#
#
# The models build in this script require the following data:
# - logAge: numeric type with log(post conception age in days, base=10). In (1),
#           we use a conversion factor of 325.25 days/year and a post conception 
#           offset of 280 days if post conception age is not available.
# - sex: factor with M or F values.
# - phenotype: numeric type with the measurements of the phenotype of interest
#
# For FreeSurfer (FS) neuroimaging phenotypes, the measurement SurfaceHoles is  
# included in the GAMLSS model. SynthSeg (SS) does not produce SurfaceHoles.
#-------------------------------------------------------------------------------

gc()
library(ggplot2)
library(gamlss) #to fit model
library(mgcv) # helps with the gam models
library(tidymv) # helps with the gam models

setwd("/Users/youngjm/Projects/mpr_analysis/r/") # change this to your path to this repo
source("lib_mpr_analysis.r")


# Load your data
fnInFs <- "/Users/youngjm/Data/slip/fs6_stats/07_fully_filtered_postcombat_clip_fs.csv"
fnInSs <- "/Users/youngjm/Data/slip/fs6_stats/07_fully_filtered_postcombat_clip_ss.csv"

dfFs <- read.csv(fnInFs)
dfSs <- read.csv(fnInSs)


# Make sure sex is a factor
dfFs$sex <- as.factor(dfFs$sex)
dfSs$sex <- as.factor(dfSs$sex)

# Make sure the dataframe has a logAge column as described at the top of the script


# Build the growth chart model
p <- "TCV" # specify the phenotype
# Incorporate SurfaceHoles in the FreeSurfer (FS) model
formulaFs <- as.formula(paste0(p, "~fp(logAge, npoly=3) + SurfaceHoles + sex - 1"))
growthChartModelFs <-gamlss(formula = formulaFs,
                            sigma.formula = formulaFs,
                            nu.formula = as.formula(paste0(p, "~1")),
                            family = GG,
                            data = na.omit(dfFs),
                            control = gamlss.control(n.cyc = 200),  # See (2)
                            trace = F)

# No SurfaceHoles in the SynthSeg (SS) model
formulaSs <- as.formula(paste0(p, "~fp(logAge, npoly=3) + sex - 1"))
growthChartModelSs <-gamlss(formula = formulaSs,
                            sigma.formula = formulaSs,
                            nu.formula = as.formula(paste0(p, "~1")),
                            family = GG,
                            data = na.omit(dfSs),
                            control = gamlss.control(n.cyc = 200),  # See (2)
                            trace = F)


# Predict the median centile of each model
medianCentileFs <- predictCentilesForAgeRange(growthChartModelFs, dfFs$logAge, 
                                              euler=median(dfFs$SurfaceHoles))
medianCentileSs <- predictCentilesForAgeRange(growthChartModelSs, dfSs$logAge)

# Calculate the age at peak (median) phenotype value
ageAtPeakFs <- 10^(sort(dfFs$logAge)[which.max(medianCentileFs)])
ageAtPeakSs <- 10^(sort(dfSs$logAge)[which.max(medianCentileSs)])


# Predict a set of centiles for each model
centileCurvesFs <- c()
centileCurvesSs <- c()

desiredCentiles <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
for (i in c(1:length(desiredCentiles))){
  centileCurvesFs[[i]] <- predictCentilesForAgeRange(growthChartModelFs, dfFs$logAge, 
                                                     euler=median(dfFs$SurfaceHoles), 
                                                     cent=desiredCentiles[[i]])
  centileCurvesSs[[i]] <- predictCentilesForAgeRange(growthChartModelSs, dfFs$logAge, 
                                                     cent=desiredCentiles[[i]])
}


# Set up a list of tick marks to use on log(post-conception age) x-axes
tickMarks <- c()
for (year in c(0, 1, 2, 5, 10, 20)){ # years
  tickMarks <- append(tickMarks, log(year*365.25 + 280, base=10))
}
tickLabels <- c("Birth", "1", "2", "5", "10", "20")


# Plot the original data and the set of centile curves on a figure
plotFs <- ggplot() +
  geom_point(aes(x=dfFs$logAge, dfFs[, p]), alpha=0.5) +
  geom_line(aes(x=dfFs$logAge, y=centileCurvesFs[[1]]), alpha=0.4) +
  geom_line(aes(x=dfFs$logAge, y=centileCurvesFs[[2]]), alpha=0.6) +
  geom_line(aes(x=dfFs$logAge, y=centileCurvesFs[[3]]), alpha=0.8) +
  geom_line(aes(x=dfFs$logAge, y=centileCurvesFs[[4]])) +
  geom_line(aes(x=dfFs$logAge, y=centileCurvesFs[[5]]), alpha=0.8) +
  geom_line(aes(x=dfFs$logAge, y=centileCurvesFs[[6]]), alpha=0.6) +
  geom_line(aes(x=dfFs$logAge, y=centileCurvesFs[[7]]), alpha=0.4) +
  scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(dfFs$logAge))) +
  labs(title=paste0("Sample Growth Chart for ", p)) + 
  xlab("Age at scan (log(years))") +
  ylab(paste0(p, " Centile")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

print(plotFs)


plotSs <- ggplot() +
  geom_point(aes(x=dfSs$logAge, dfSs[, p]), alpha=0.5) +
  geom_line(aes(x=dfSs$logAge, y=centileCurvesSs[[1]]), alpha=0.4) +
  geom_line(aes(x=dfSs$logAge, y=centileCurvesSs[[2]]), alpha=0.6) +
  geom_line(aes(x=dfSs$logAge, y=centileCurvesSs[[3]]), alpha=0.8) +
  geom_line(aes(x=dfSs$logAge, y=centileCurvesSs[[4]])) +
  geom_line(aes(x=dfSs$logAge, y=centileCurvesSs[[5]]), alpha=0.8) +
  geom_line(aes(x=dfSs$logAge, y=centileCurvesSs[[6]]), alpha=0.6) +
  geom_line(aes(x=dfSs$logAge, y=centileCurvesSs[[7]]), alpha=0.4) +
  scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(dfSs$logAge))) +
  labs(title=paste0("Sample Growth Chart for ", p)) + 
  xlab("Age at scan (log(years))") +
  ylab(paste0(p, " Centile")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

print(plotFs)


# Calculate the closest centile to each subject's phenotype value
phenoCentilesFs <- calculatePhenotypeCentile(growthChartModelFs, dfFs[[p]], dfFs$logAge, dfFs$sex, dfFs$SurfaceHoles)
regionsFs <- rep(p, length(phenoCentilesFs))
idxesFs <- c(1:length(phenoCentilesFs))

phenoCentilesSs <- calculatePhenotypeCentile(growthChartModelSs, dfSs[[p]], dfSs$logAge, dfSs$sex)
regionsSs <- rep(p, length(phenoCentilesSs))
idxesSs <- c(1:length(phenoCentilesSs))


# Plot the centiles of the phenotypes in a violin plot
dfFsViolin <- data.frame(idxesFs, regionsFs, phenoCentilesFs)
plotPhenoCentFs <- ggplot(data=dfFsViolin, aes(regionsFs, phenoCentilesFs)) +
  geom_violin(color="gray", fill="gray", alpha=0.35) +
  geom_jitter(height = 0, width=0.15, alpha=0.65) +
  labs(title=paste0("Centiles (Phenotype = ", p,")")) +
  xlab("FreeSurfer") +
  ylab("Centile Value") +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

print(plotPhenoCentFs)


dfSsViolin <- data.frame(idxesSs, regionsSs, phenoCentilesSs)
plotPhenoCentSs <- ggplot(data=dfSsViolin, aes(regionsSs, phenoCentilesSs)) +
  geom_violin(color="gray", fill="gray", alpha=0.35) +
  geom_jitter(height = 0, width=0.15, alpha=0.65) +
  labs(title=paste0("Centiles (Phenotype = ", p,")")) +
  xlab("SynthSeg") +
  ylab("Centile Value") +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

print(plotPhenoCentSs)
