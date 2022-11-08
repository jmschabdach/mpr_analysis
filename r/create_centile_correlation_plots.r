gc()
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork) # graph organization within a figure

library(gamlss) #to fit model
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")


fn <- '/Users/youngjm/Data/lifespan_growth_charts/lifespan_centile_medians.csv'
lc <- read.csv(fn)
fn <- '/Users/youngjm/Data/clip/fs6_stats/06_combatted_plus_metadata.csv'
clipDf <- read.csv(fn)
fn <- '/Users/youngjm/Data/clip/fs6_stats/synthseg_2.0_phenotypes_cleaned.csv'
ssDf <- read.csv(fn) 

names(ssDf)[names(ssDf) == "subj_id"] <- 'patient_id'

# Prep Lifespan centiles df to match CLIP formatting for sex
lc <- lc %>%
  mutate(sex = str_replace(sex, 'Female', 'F'))
lc <- lc %>%
  mutate(sex = str_replace(sex, 'Male', 'M'))

# Add a column to Lifespan for log(age)
lc$logAge <- log(lc$age, 10)

# Scale the values - they're N x 10,000
lc$value <- lc$value * 10000

# Add a column to synthseg for post conception age
ssDf$age <- ssDf$age_at_scan_days+280
ssDf$logAge <- log(ssDf$age, 10)
ssDf <- addPrimaryScanReasonCol(ssDf)
ssDf$scanner_id <- as.factor(ssDf$scanner_id)
ssDf$sex <- as.factor(ssDf$sex)

# Rename the columns of the combatted clip dataframe we care about
names(ssDf)[names(ssDf) == "Cortex"] <- 'GMV'
names(ssDf)[names(ssDf) == "Ventricles"] <- 'CSF'

# Replace any 0 valued phenotypes with 0.001
ssDf[, "WMV"][ssDf[, "WMV"] == 0] <- NA
ssDf[, "CSF"][ssDf[, "CSF"] == 0] <- NA


# Add a column to clip for post conception age
clipDf$age <- clipDf$age_at_scan_days+280
clipDf <- addPrimaryScanReasonCol(clipDf)

# Add a column to clip for log(age)
clipDf$logAge <- log(clipDf$age, 10)

# Rename the columns of the combatted clip dataframe we care about
names(clipDf)[names(clipDf) == "TotalGrayVol"] <- 'GMV'
names(clipDf)[names(clipDf) == "CerebralWhiteMatterVol"] <- 'WMV'
names(clipDf)[names(clipDf) == "SubCortGrayVol"] <- 'sGMV'
names(clipDf)[names(clipDf) == "VentricleVolume"] <- 'CSF'
names(clipDf)[names(clipDf) == "CorticalSurfaceArea"] <- 'SA'
names(clipDf)[names(clipDf) == "MeanCorticalThickness"] <- 'CT'
# names(clipDf)[names(clipDf) == "TCVTransformed.normalised"] <- 'TCV'
clipDf$sex <- as.factor(clipDf$sex)
# clipDf <- clipDf[clipDf$age < max(clipDf$age), ]
clipDf$scanner_id <- as.factor(clipDf$scanner_id)

# Filter out any subjects not in both FS and SS dataframes
ssDf <- drop_na(ssDf)
ssDf <- ssDf[(ssDf$patient_id %in% clipDf$patient_id), ]
clipDf <- clipDf[(clipDf$patient_id %in% ssDf$patient_id), ]

phenos <- c('GMV', 'WMV', 'sGMV', 'CSF', 'TCV') # 'CT', 'SA', 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

minAgeDaysPostCon <- min(clipDf$logAge)
maxAgeDaysPostCon <- max(clipDf$logAge)
ageLimited <- unique(lc[lc$logAge > minAgeDaysPostCon & lc$logAge < maxAgeDaysPostCon, "logAge"])

r <- c()
rSs <- c()
rFsSs <- c()
rFsSsCent <- c()
ageAtPeakCLIP <- c()
ageAtPeakLifespan <- c()
ageAtPeakSs <- c()
centiles <- c()
regions <- c()
idxes <- c()
reasons <- c()
bigPlots <- c()

# Calculating the centile for a subject
calculatePhenotypeCentile <- function(model, measuredPhenotypeValue, logAge, surfaceHoles, sex){
  centileDistribution <- 1:9999/10000
  centiles <- c()
  for (i in 1:length(measuredPhenotypeValue)){
    newData <- data.frame(logAge=logAge[[i]],
                          SurfaceHoles=surfaceHoles[[i]],
                          sex=sex[[i]])
    predModel <- predictAll(model, newdata=newData)
    expectedPhenotypeValue <- qGG(centileDistribution, mu=predModel$mu, sigma=predModel$sigma, nu=predModel$nu)
    centiles[i] <- centileDistribution[which.min(abs(measuredPhenotypeValue[[i]] - expectedPhenotypeValue))]
  }
  return(centiles)
}

calculatePhenotypeCentileSynthSeg <- function(model, measuredPhenotypeValue, logAge, surfaceHoles, sex){
  centileDistribution <- 1:9999/10000
  centiles <- c()
  for (i in 1:length(measuredPhenotypeValue)){
    newData <- data.frame(logAge=logAge[[i]],
                          sex=sex[[i]])
    predModel <- predictAll(model, newdata=newData)
    expectedPhenotypeValue <- qGG(centileDistribution, mu=predModel$mu, sigma=predModel$sigma, nu=predModel$nu)
    centiles[i] <- centileDistribution[which.min(abs(measuredPhenotypeValue[[i]] - expectedPhenotypeValue))]
  }
  return(centiles)
}

for ( p in phenos ) {
  print(p)
  plots <- c()
  fnOut <- paste0('/Users/youngjm/Data/clip/figures/2022-10-19_clip_lifespan_correlation_', p, '.png')
  smallLc <- lc[lc$feature==p, ]
  
  # For plot 2:
  # 1. Generate GAMLSS models
  formula <- as.formula(paste0(p, "~fp(logAge, npoly=3) + SurfaceHoles + sex - 1"))
  gamModel <-gamlss(formula = formula,
                    sigma.formula = formula,
                    nu.formula = as.formula(paste0(p, "~1")),
                    family = GG,
                    data = na.omit(clipDf),
                    # data = rbind(na.omit(clipDf), placeholderDf),
                    control = gamlss.control(n.cyc = 200),  # lifespan
                    trace = F)
  print("finished training the models")

  # 2. Predict phenotype values for set age range
  newDataM <- data.frame(logAge=sort(ageLimited),
                         SurfaceHoles=c(rep(median(clipDf$SurfaceHoles), length(ageLimited))),
                         sex=c(rep(as.factor("M"),  length(ageLimited))))
  clipPredModelM <- predictAll(gamModel, newdata=newDataM)
  
  newDataF <- data.frame(logAge=sort(ageLimited),
                         SurfaceHoles=c(rep(median(clipDf$SurfaceHoles), length(ageLimited))),
                         sex=c(rep(as.factor("F"),  length(ageLimited))))
  clipPredModelF <- predictAll(gamModel, newdata=newDataF)

  # The c(0.5) is for the 50th percentile
  phenoMedianPredsM <- qGG(c(0.5), 
                          mu=clipPredModelM$mu, 
                          sigma=clipPredModelM$sigma, 
                          nu=clipPredModelM$nu)

  phenoMedianPredsF <- qGG(c(0.5), 
                           mu=clipPredModelF$mu, 
                           sigma=clipPredModelF$sigma, 
                           nu=clipPredModelF$nu)
  phenoMedianPreds <- (phenoMedianPredsF + phenoMedianPredsM)/2

  # 3. Get the smaller Lifespan data frame
  smallLc <- lc[(lc$logAge %in% ageLimited) & (lc$feature == p) & (lc$sex == 'M'), ]
  smallLc$sex <- as.factor(smallLc$sex)
  
  
  # 4. Calculate correlation, age at peak, and predicted centiles
  r[[p]] <- cor(smallLc$value, phenoMedianPreds)
  ageAtPeakCLIP[[p]] <- 10^(sort(ageLimited)[which.max(phenoMedianPreds)])
  ageAtPeakLifespan[[p]] <- smallLc[smallLc$value == max(smallLc$value), 'age']
  tmp <- calculatePhenotypeCentile(gamModel, clipDf[[p]], clipDf$logAge, clipDf$SurfaceHoles, clipDf$sex)
  centiles <- append(centiles, tmp)
  regions <- append(regions, rep(p, length(tmp)))
  idxes <- append(idxes, c(1:length(tmp)))
  reasons <- append(reasons, clipDf$top_scan_reason_factors)
  
  # 5. Plot CLIP vs Lifespan
  plots[[1]] <- ggplot() +
    geom_point(data=clipDf, aes(x=logAge, y=clipDf[,p], color=top_scan_reason_factors), alpha=0.3) +
    geom_line(aes(x=ageLimited, y=phenoMedianPreds, linetype="Predicted for CLIP")) +
    # geom_line(aes(x=ageLimited, y=phenoMedianPredsM, linetype="Predicted for CLIP M"), color="blue") +
    # geom_line(aes(x=ageLimited, y=phenoMedianPredsF, linetype="Predicted for CLIP F"), color="orange") +
    geom_line(data=smallLc, aes(x=logAge, y=value, linetype="Lifespan")) +
    # geom_text(x=(0.85*max(clipDf$logAge)), y=(0.85*max(clipDf[,p])), label=paste0("r=", format(r[[p]], digits=4), "(50th Centiles)")) +
    scale_color_manual(values = cbbPalette, name = "Reason for Scan") +
    scale_linetype_manual(values = c('solid', 'dashed', 'solid', 'solid'), name="50th Centile")+
    theme(plot.title=element_text(hjust=0.5)) +
    labs(subtitle = paste0("FreeSurfer v. 6.0.0 (r=", format(r[[p]], digits=4), ")")) +
    xlab("log(Age) (days)") +
    # xlim(0, max(clipDf$logAge)) +
    # ylim(0, max(clipDf[,p])) +
    ylab("Phenotype Value") + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  w <- 320
  
  if (p %in% colnames(ssDf)){
    # Slightly different model for SynthSeg
    formulaSs <- as.formula(paste0(p, "~fp(logAge, npoly=3) + sex - 1"))
    gamModelSs <-gamlss(formula = formulaSs,
                        sigma.formula = formulaSs,
                        nu.formula = as.formula(paste0(p, "~1")),
                        family = GG,
                        data = na.omit(ssDf),
                        # data = rbind(na.omit(clipDf), placeholderDf),
                        control = gamlss.control(n.cyc = 200),  # lifespan
                        trace = F)
    
    # 2. Predict phenotype values for set age range
    newDataM <- data.frame(logAge=sort(ageLimited),
                           sex=c(rep(as.factor("M"),  length(ageLimited))))
    clipPredModelSsM <- predictAll(gamModelSs, newdata=newDataM)
    
    newDataF <- data.frame(logAge=sort(ageLimited),
                           sex=c(rep(as.factor("F"),  length(ageLimited))))
    clipPredModelSsF <- predictAll(gamModelSs, newdata=newDataF)
    
    # The c(0.5) is for the 50th percentile
    phenoMedianPredsMSs <- qGG(c(0.5), 
                             mu=clipPredModelSsM$mu, 
                             sigma=clipPredModelSsM$sigma, 
                             nu=clipPredModelSsM$nu)
    
    phenoMedianPredsFSs <- qGG(c(0.5), 
                             mu=clipPredModelSsF$mu, 
                             sigma=clipPredModelSsF$sigma, 
                             nu=clipPredModelSsF$nu)
    phenoMedianPredsSs <- (phenoMedianPredsFSs + phenoMedianPredsMSs)/2

    rSs[[p]] <- cor(smallLc$value, phenoMedianPredsSs)
    rFsSs[[p]] <- cor(phenoMedianPreds, phenoMedianPredsSs)
    rFsSsCent[[p]] <- cor(clipDf[,p], ssDf[,p])
    ageAtPeakSs[[p]] <- 10^(sort(ageLimited)[which.max(phenoMedianPredsSs)])

    plots[[2]] <- ggplot(x="linear") +
      geom_point(data=ssDf, aes(x=logAge, y=ssDf[,p], color=top_scan_reason_factors), alpha=0.3) +
      geom_line(aes(x=ageLimited, y=phenoMedianPredsSs, linetype="Predicted for CLIP")) +
      # geom_line(aes(x=ageLimited, y=phenoMedianPredsMSs, linetype="Predicted for CLIP M"), color="blue") +
      # geom_line(aes(x=ageLimited, y=phenoMedianPredsFSs, linetype="Predicted for CLIP F"), color="orange") +
      geom_line(data=smallLc, aes(x=logAge, y=value, linetype="Lifespan")) +
      scale_color_manual(values = cbbPalette, name = "Reason for Scan") +
      scale_linetype_manual(values = c('solid', 'dashed', 'solid', 'solid'), name="50th Centile")+
      labs(subtitle = paste0("SynthSeg 2.0 (r=", format(rSs[[p]], digits=4), ")")) +
      xlab("log(Age) (days)") +
      # xlim(0, max(ssDf$logAge)) +
      # ylim(0, max(ssDf[,p])) +
      ylab("Phenotype Value") + 
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
    
    plots[[3]] <- ggplot() +
      geom_point(aes(x=clipDf[,p], y=ssDf[,p], color=clipDf$top_scan_reason_factors), alpha=0.3) +
      # geom_abline(slope = rFsSs[[p]]) +
      # geom_smooth(aes(x=clipDf[,p], y=ssDf[,p]), color="red") +
      scale_color_manual(values = cbbPalette, name = "Reason for Scan") +
      labs(subtitle = paste0(p, ": FreeSurfer vs SynthSeg (r=", format(rFsSs[[p]], digits=4), ")")) +
      xlab("FreeSurfer 6.0.0") +
      ylab("SynthSeg 2.0") + 
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

    w <- 900
  }
  
  # Show the plots
  patch <- wrap_plots(plots, guides='collect')
  # bigPlots <- append(bigPlots, plots)
  png(file=fnOut,
      width=w, height=250)
  print(patch + plot_annotation(title=paste0("Comparison of Lifespan and CLIP for ", p)))
  dev.off()
}

# Make a regular centile plot
smallLc <- lc[lc$feature=="GMV", ]

# For plot 2:
# 1. Generate GAMLSS models
formula <- as.formula("GMV~fp(logAge, npoly=3) + SurfaceHoles + sex - 1")
gamModel <-gamlss(formula = formula,
                  sigma.formula = formula,
                  nu.formula = as.formula("GMV~1"),
                  family = GG,
                  data = na.omit(clipDf),
                  # data = rbind(na.omit(clipDf), placeholderDf),
                  control = gamlss.control(n.cyc = 200),  # lifespan
                  trace = F)
print("finished training the models")

# 2. Predict phenotype values for set age range
newDataM <- data.frame(logAge=sort(ageLimited),
                       SurfaceHoles=c(rep(median(clipDf$SurfaceHoles), length(ageLimited))),
                       sex=c(rep(as.factor("M"),  length(ageLimited))))
clipPredModelM <- predictAll(gamModel, newdata=newDataM)

newDataF <- data.frame(logAge=sort(ageLimited),
                       SurfaceHoles=c(rep(median(clipDf$SurfaceHoles), length(ageLimited))),
                       sex=c(rep(as.factor("F"),  length(ageLimited))))
clipPredModelF <- predictAll(gamModel, newdata=newDataF)

# The c(0.5) is for the 50th percentile
fanCentiles <- c()
desiredCentiles <- c(0.004, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.996)
for (i in c(1:length(desiredCentiles))){
  print(desiredCentiles[[i]])
  print(i)
  phenoMedianPredsM <- qGG(desiredCentiles[[i]], 
                           mu=clipPredModelM$mu, 
                           sigma=clipPredModelM$sigma, 
                           nu=clipPredModelM$nu)
  
  phenoMedianPredsF <- qGG(desiredCentiles[[i]], 
                           mu=clipPredModelF$mu, 
                           sigma=clipPredModelF$sigma, 
                           nu=clipPredModelF$nu)
  fanCentiles[[i]] <- (phenoMedianPredsF + phenoMedianPredsM)/2
}

sampleCentileFan <- ggplot() +
  geom_point(aes(x=clipDf$logAge, clipDf$GMV)) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[1]]), alpha=0.2) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[2]]), alpha=0.4) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[3]]), alpha=0.6) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[4]]), alpha=0.8) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[5]])) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[6]]), alpha=0.8) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[7]]), alpha=0.6) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[8]]), alpha=0.4) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[9]]), alpha=0.2) +
  labs(title="Sample Centile Growth Chart for GMV") + 
  xlab("Age at Scan") +
  ylab("GMV") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Make the violin plots
violinDf <- data.frame(idxes, regions, centiles, reasons)
violin <- ggplot(data=violinDf, aes(regions, centiles)) +
  geom_violin(color="gray", fill="gray", alpha=0.5) +
  geom_jitter(height = 0, width=0.15, aes(color=reasons), alpha=0.35) +
  scale_color_manual(values = cbbPalette, name = "Reason for Scan") +
  labs(title="Distributions of Predicted Centiles (FS)") + 
  xlab("Tissue Type") +
  ylab("Predicted Centile") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

reasonPlots <- c()
uniqueReasons <- levels(violinDf$reasons)
for (r in c(1:length(uniqueReasons))){
  reasonDf <- violinDf[reasons==uniqueReasons[[r]], ]
  reasonPlots[[r]] <- ggplot(data=reasonDf, aes(regions, centiles)) +
    geom_violin(color="gray", fill="gray", alpha=0.35) +
    geom_jitter(height = 0, width=0.15, aes(color=reasons), alpha=0.45) +
    scale_color_manual(values = cbbPalette[[r]]) +
    labs(title=paste0("Centiles (Reason = ", uniqueReasons[[r]],")")) + 
    xlab("Tissue Type") +
    ylab("Predicted Centile") + 
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}
topPanel <- wrap_plots(sampleCentileFan + violin)
bottomPanel <- wrap_plots(reasonPlots, ncol=3, guides="collect")
layout <- layout <-"
A
B
B
"
patch <- wrap_plots(topPanel + bottomPanel + plot_layout(design=layout), guides="auto")
png(file="/Users/youngjm/Data/clip/figures/2022-11-08_clip_predicted_centiles.png",
    width=1100, height=800)
print(patch)
dev.off()

convertToYears <- function(ages){
  for (i in names(ages)){
    print(i)
    print((ages[[i]]-280)/365.25)
  }
}

convertToYears(ageAtPeakCLIP)
convertToYears(ageAtPeakSs)
convertToYears(ageAtPeakLifespan)

# Statistical tests: for each phenotype, compare the centile values across reason for scan (ANOVA)
# centile ~ reason for scan | phenotype
violinDf$reasons <- as.factor(violinDf$reasons)
for (p in phenos){
  aovOut <- aov(centiles ~ reasons, data=na.omit(violinDf[regions == p,]))
  print(summary(aovOut))
}
