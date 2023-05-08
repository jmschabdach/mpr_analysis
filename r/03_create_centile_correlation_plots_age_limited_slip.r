gc()
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)
library(stringr)
library(patchwork) # graph organization within a figure

library(gamlss) #to fit model
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")

# Part 0: load the data from multiple files and make sure all of the data
#  types match
fn <- '/Users/youngjm/Data/lifespan_growth_charts/slip-agemedian_centiles.csv'
lc <- read.csv(fn)
fn <- '/Users/youngjm/Data/slip/fs6_stats/06_combatted_fs_plus_metadata.csv'
clipDf <- read.csv(fn)
fn <- '/Users/youngjm/Data/slip/fs6_stats/06_combatted_ss_plus_metadata.csv'
ssDf <- read.csv(fn) 
preCombatFn <- '/Users/youngjm/Data/slip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
preCombatDf <- read.csv(preCombatFn)

# SLIP age lc
lc$feature <- as.factor(lc$feature)
names(lc)[names(lc) == "age"] <- "logAge"
lc$age <- (10^lc$logAge) # post conception age in days

# Prep synthseg data
# Add a column to synthseg for post-conception age
ssDf$age <- ssDf$age_at_scan_days+280
# Add log of post-conception age
ssDf$logAge <- log(ssDf$age, 10)
# Convert scanner id to a factor
ssDf$scanner_id <- gsub("Scanner ","",as.character(ssDf$scanner_id))
ssDf$scanner_id <- as.factor(ssDf$scanner_id)
ssDf$sex <- as.factor(ssDf$sex)
ssDf$proc_ord_year <- as.factor(ssDf$proc_ord_year)
ssDf$top_scan_reason_factors <- as.factor(ssDf$top_scan_reason_factors)
# Rename columns to match the clipDf columns
names(ssDf)[names(ssDf) == "Cortex"] <- 'GMV'
names(ssDf)[names(ssDf) == "Ventricles"] <- 'CSF'
names(ssDf)[names(ssDf) == "subj_id"] <- 'patient_id'
# Replace any 0 valued phenotypes with NA
ssDf[, "WMV"][ssDf[, "WMV"] <= 0] <- NA
ssDf[, "CSF"][ssDf[, "CSF"] <= 0] <- NA

# Prep CLIP FS data
# Add a column to clip for log(post conception age)
clipDf$logAge <- log(clipDf$age_at_scan_days+280, 10)
# Make sex and scanner_id factors
clipDf$sex <- as.factor(clipDf$sex)
clipDf$scanner_id <-gsub("Scanner ","",as.character(clipDf$scanner_id))
clipDf$scanner_id <- as.factor(clipDf$scanner_id)
clipDf$proc_ord_year <- as.factor(clipDf$proc_ord_year)
clipDf$top_scan_reason_factors <- as.factor(clipDf$top_scan_reason_factors)
# Drop rows with missing values
clipDf <- drop_na(clipDf)

# Prep the precombat FS data
preCombatDf$sex <- as.factor(preCombatDf$sex)
preCombatDf$logAge <- log((preCombatDf$age_at_scan_days + 280), base=10)
preCombatDf$scanner_id <- gsub("Scanner ","",as.character(preCombatDf$scanner_id))
preCombatDf$scanner_id <- as.factor(preCombatDf$scanner_id)
preCombatDf$top_scan_reason_factors <- as.factor(preCombatDf$top_scan_reason_factors)
preCombatDf$CSF <- abs(preCombatDf$BrainSegVol - preCombatDf$BrainSegVolNotVent)
preCombatDf <- drop_na(preCombatDf)

# Do an inner join on the pre and post combat FS dataframes
clipDf <- inner_join(clipDf, preCombatDf, 
           by=c("patient_id", "scanner_id", "age_at_scan_days", "logAge", "top_scan_reason_factors"), 
           suffix=c("", ".pre"))

# Filter out any subjects not in both FS and SS dataframes
ssDf <- drop_na(ssDf)
ssDf <- ssDf[(ssDf$patient_id %in% clipDf$patient_id), ]
clipDf <- clipDf[(clipDf$patient_id %in% ssDf$patient_id), ]

# Adding another area for QC: abs(x-y)/max(x,y)
setorder(clipDf, age_at_scan_days)  
setorder(ssDf, age_at_scan_days)

# Set up variables needed for the analysis
phenos <- c('GMV', 'WMV', 'sGMV', 'CSF', 'TCV') # 'CT', 'SA', 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",  
                "#D55E00", "#CC79A7", "#0072B2", "#F0E442")

# Determine the limited age range for predicting the centiles with the GAMLSS
minAgeDaysPostCon <- min(clipDf$logAge)
maxAgeDaysPostCon <- max(clipDf$logAge)
ageLimited <- unique(lc[lc$logAge > minAgeDaysPostCon & lc$logAge < maxAgeDaysPostCon, "logAge"])
lc <- lc[lc$logAge > minAgeDaysPostCon & lc$logAge < maxAgeDaysPostCon, ]

write.csv(clipDf, "/Users/youngjm/Data/slip/fs6_stats/07_fully_filtered_postcombat_clip_fs.csv")

# Set up lists for storing variables in the following for loops
r <- c()
rSs <- c()
rFsSs <- c()
rFsSsCent <- c()
ageAtPeakCLIP <- c()
ageAtPeakLifespan <- c()
ageAtPeakSs <- c()
centiles <- c()
preCombatCentiles <- c()
phenoDist <- c()
preCombatPhenoDist <- c()
regions <- c()
idxes <- c()
reasons <- c()
yearOfScan <- c()
bigPlots <- c()

reasonLabels <- c("Developmental disorder", "Clinical eye or vision finding", "Headache", "Other", "Suspected seizure")

# Set up a list of tick marks to use on log(post-conception age) x-axes
tickMarks <- c()
for (year in c(0, 1, 2, 5, 10, 20)){ # years
  tickMarks <- append(tickMarks, log(year*365.25 + 280, base=10))
}
tickLabels <- c("Birth", "1", "2", "5", "10", "20")


# Generate GAMLSS models for each phenotype
for ( p in phenos ) {
  print(p)
  plots <- c() # reset the list of plots
  fnOut <- paste0('/Users/youngjm/Data/slip/figures/2022-10-19_clip_agelimited_lifespan_correlation_', p, '.png')

  # Pull the pre and post combat phenotypes
  phenoDist[[p]] <- clipDf[,p]
  preCombatPhenoDist[[p]] <- clipDf[,paste0(p, ".pre")]

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

  # Predict the median centile for the M and F models
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
  smallLc <- lc[lc$feature==p, ]
  
  
  
  # 4. Calculate correlation, age at peak, and predicted centiles
  r[[p]] <- cor(smallLc$value, phenoMedianPreds)
  ageAtPeakCLIP[[p]] <- 10^(sort(ageLimited)[which.max(phenoMedianPreds)])
  ageAtPeakLifespan[[p]] <- smallLc[smallLc$value == max(smallLc$value), 'age']
  tmp <- calculatePhenotypeCentile(gamModel, clipDf[[p]], clipDf$logAge, clipDf$SurfaceHoles, clipDf$sex)
  centiles[[p]] <- tmp
  preCombatCentiles[[p]] <- calculatePhenotypeCentile(gamModel, clipDf[[paste0(p, ".pre")]], clipDf$logAge, clipDf$SurfaceHoles, clipDf$sex)
  regions <- append(regions, rep(p, length(tmp)))
  idxes <- append(idxes, c(1:length(tmp)))
  reasons <- append(reasons, clipDf$top_scan_reason_factors)
  yearOfScan <- append(yearOfScan, clipDf$proc_ord_year)
  
  
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
  
  # Get the y min and y max values
  ymin <- min(min(clipDf[,p]), min(ssDf[,p]))
  ymax <- max(max(clipDf[,p]), max(ssDf[,p]))
  
  # 5. Plot CLIP vs Lifespan
  plots[[1]] <- ggplot() +
    geom_point(data=clipDf, aes(x=logAge, y=clipDf[,p], color=top_scan_reason_factors), alpha=0.3) +
    geom_line(aes(x=ageLimited, y=phenoMedianPreds, linetype="Predicted for SLIP")) +
    geom_line(data=smallLc, aes(x=logAge, y=value, linetype="SLIP-age LBCC")) +
    scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
    scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'solid'), name="50th Centile")+
    scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                       limits=c(tickMarks[[1]], max(clipDf$logAge))) +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(subtitle = paste0("FreeSurfer v. 6.0.0 (r=", format(r[[p]], digits=3), ")")) +
    xlab("Age at scan (log(years))") +
    ylab("Phenotype Value") + 
    ylim(ymin, ymax) +
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 18))

  plots[[2]] <- ggplot(x="linear") +
    geom_point(data=ssDf, aes(x=logAge, y=ssDf[,p], color=top_scan_reason_factors), alpha=0.3) +
    geom_line(aes(x=ageLimited, y=phenoMedianPredsSs, linetype="Predicted for SLIP")) +
    geom_line(data=smallLc, aes(x=logAge, y=value, linetype="SLIP-age LBCC")) +
    scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
    scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'solid'), name="50th Centile")+
    labs(subtitle = paste0("SynthSeg+ (r=", format(rSs[[p]], digits=3), ")")) +
    xlab("Age at scan (log(years))") +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                       limits=c(tickMarks[[1]], max(clipDf$logAge))) +
    ylab("Phenotype Value") + 
    ylim(ymin, ymax) +
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 18))
  
  plots[[3]] <- ggplot() +
    geom_point(aes(x=clipDf[,p], y=ssDf[,p], color=clipDf$top_scan_reason_factors), alpha=0.3) +
    scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
    labs(subtitle = paste0("FreeSurfer vs SynthSeg (r=", format(rFsSsCent[[p]], digits=3), ")")) +
    xlab("FreeSurfer 6.0.0") +
    ylab("SynthSeg+") + 
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 18))

  w <- 1400
  
  # Show the plots
  patch <- wrap_plots(plots, nrow = 1, guides='collect')
  png(file=fnOut,
      width=w, height=300)
  print(patch + plot_annotation(title=p, theme = theme(text=element_text(size=18))))
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
desiredCentiles <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)

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
  geom_point(aes(x=clipDf$logAge, clipDf$GMV), alpha=0.5) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[1]]), alpha=0.4) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[2]]), alpha=0.6) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[3]]), alpha=0.8) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[4]])) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[5]]), alpha=0.8) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[6]]), alpha=0.6) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[7]]), alpha=0.4) +
  scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(clipDf$logAge))) +
  labs(title="Sample Centile Growth Chart for GMV") + 
  xlab("Age at scan (log(years))") +
  ylab("GMV Centile") + 
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

rawPhenoPlot <- ggplot(x="linear") +
  geom_point(data=clipDf, aes(x=logAge, y=GMV, color=top_scan_reason_factors), alpha=0.3) +
  scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
  labs(title="Phenotype Distribution for GMV") + 
  xlab("Age at scan (log(years))") +
  ylab("GMV Value") + 
  scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(clipDf$logAge))) +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 18))


# Make the violin plots
centilesList <- c(centiles$GMV, centiles$WMV, centiles$sGMV, centiles$CSF, centiles$TCV)
vScanId <- c(rep(clipDf$scan_id, 5))
violinDf <- data.frame(idxes, regions, centilesList, reasons, yearOfScan, vScanId)
write.csv(violinDf, "/Users/youngjm/Data/slip/fs6_stats/fs_gamlss_centiles.csv", row.names = FALSE)
violin <- ggplot(data=violinDf, aes(regions, centilesList)) +
  geom_violin(color="gray", fill="gray", alpha=0.5) +
  geom_jitter(height = 0, width=0.15, aes(color=reasons), alpha=0.65) +
  scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
  labs(title="Distributions of Centiles (FS)") + 
  xlab("Tissue Type") +
  ylab("Centile") + 
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 18))

# Change this section from scan reason focused to phenotype focused
phenoCentilePlots <- c()
phenoCentilePlots[[1]] <- violin
for (p in c(1:length(phenos))){
  phenoDf <- violinDf[regions==phenos[[p]], ]
  phenoCentilePlots[[p+1]] <- ggplot(data=phenoDf, aes(reasons, centilesList)) +
    geom_violin(color="gray", fill="gray", alpha=0.35) +
    geom_jitter(height = 0, width=0.15, aes(color=reasons), alpha=0.65) +
    scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
    labs(title=paste0("Centile (Phenotype = ", phenos[[p]],")")) + 
    xlab("Reason for Scan") +
    scale_x_discrete(labels=c("DD", "EV", "H", "O", "SS")) +
    ylab("Centile") + 
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 20))
}
topPanel <- wrap_plots(rawPhenoPlot + sampleCentileFan)
bottomPanel <- wrap_plots(phenoCentilePlots, ncol=2, guides="collect")
layout <-"
A
B
B
B
"
patch <- wrap_plots(topPanel + bottomPanel + plot_layout(design=layout), guides="auto")
png(file="/Users/youngjm/Data/slip/figures/2022-11-08_clip_predicted_centiles.png",
    width=1400, height=1400)
print(patch)
dev.off()



convertAgeToYears(ageAtPeakCLIP)
convertAgeToYears(ageAtPeakSs)
convertAgeToYears(ageAtPeakLifespan)

# Statistical tests: for each phenotype, compare the centile values across reason for scan (ANOVA)
# centile ~ reason for scan | phenotype
violinDf$reasons <- as.factor(violinDf$reasons)
violinDf$yearOfScan <- as.factor(yearOfScan)
for (p in phenos){
  aovOut <- aov(centilesList ~ yearOfScan, data=na.omit(violinDf[regions == p,]))
  # lm , yearOfScan = numeric
    print(summary(aovOut))
}



# Test the distributions of phenotype values before and after ComBat
plots <- c()
pvals <- c()
for (i in seq(1:length(phenos))){
  print(i)
  p <- phenos[[i]]
  # Set up the data frame
  tmp <- data.frame(PreCombat= preCombatCentiles[[i]],
                    PostCombat = centiles[[i]],
                    scannerId = clipDf$scanner_id)
  reps <- length(tmp$PreCombat)

  # Transform the data
  tmp <- data.frame(Time = c(rep("pre-ComBat", reps), rep("post-ComBat", reps)),
                    Phenotype=c(tmp$PreCombat, tmp$PostCombat),
                    Scanner = c(tmp$scannerId, tmp$scannerId))

  # Reorder Time
  tmp$Time <- factor(tmp$Time, levels = c("pre-ComBat", "post-ComBat"))
  tmp$Scanner <- factor(tmp$Scanner)
  for (time in levels(tmp$Time)){
    # Visualize - separate by scanner_id
    print(paste(p, time))
    miniDf <- na.omit(tmp[tmp$Time == time,])
    # Perform the test
    test <- aov(miniDf$Phenotype ~ miniDf$Scanner, data=na.omit(miniDf))
    pval <- summary(test)[[1]][["Pr(>F)"]][1]
    fval <- summary(test)[[1]][["F value"]][1]
    dof <- summary(test)[[1]][["Df"]][1]
    pvals[[paste0(p, "_", time)]] <- p
    plots[[paste0(p, "_", time)]] <- ggplot(tmp[tmp$Time == time, ]) +
      aes(x=Scanner, y=Phenotype) +
      geom_boxplot(aes(fill=Scanner), varwidth = TRUE, alpha=0.75) +
      geom_jitter(height = 0, width=0.15, alpha=0.4) +
      scale_fill_viridis_d(name = "Scanner ID") +
      labs(title=paste0(p, " ", time),
           subtitle=paste0("(F(",paste(dim(clipDf)[1]-1), ", ", paste(dof),") = ",
                        format(fval, digits=3), ", p < ", format(pval, digits=3), ")")) +
      theme(legend.position = "none",
            axis.line = element_line(colour = "black"),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 18))
  }
}

imgOut <- "~/Data/slip/figures/2022-11-16_pre_post_combat_combo_all.png"
patch <- wrap_plots(plots, ncol=2, byrow=TRUE, guides="collect")
png(file=imgOut,
    width=1000, height=1400)
print(patch)
dev.off()

# Pull the stats for the table
print(table(clipDf$sex))
clipDf$age_group <- case_when(
  clipDf$age_in_years < 2.0 ~ "0-2",
  clipDf$age_in_years < 5.0 ~ "02-5",
  clipDf$age_in_years < 10.0 ~ "05-10",
  clipDf$age_in_years < 13.0 ~ "10-13",
  clipDf$age_in_years < 18.0 ~ "13-18",
  TRUE ~ "18+"
)
clipDf$age_group <- as.factor(clipDf$age_group)

clipDf$scan_year_group <- case_when(
  clipDf$proc_ord_year %in% c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020") ~ "2011+",
  TRUE ~ as.character(clipDf$proc_ord_year)
)
clipDf$scan_year_group <- as.factor(clipDf$scan_year_group)

print(table(clipDf$age_group))
print(table(clipDf[clipDf$sex == "M", "age_group"]))
print(table(clipDf[clipDf$sex == "F", "age_group"]))

print(table(clipDf$top_scan_reason_factors))
print(table(clipDf[clipDf$sex == "M", "top_scan_reason_factors"]))
print(table(clipDf[clipDf$sex == "F", "top_scan_reason_factors"]))

print(table(clipDf$scan_year_group))
print(table(clipDf[clipDf$sex == "M", "scan_year_group"]))
print(table(clipDf[clipDf$sex == "F", "scan_year_group"]))

scanYearGroup <- case_when(
  yearOfScan %in% c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020") ~ "2011+",
  TRUE ~ as.character(yearOfScan)
)

violinDf <- data.frame(idxes, regions, centilesList, reasons, scanYearGroup)

# Change this section from scan reason focused to phenotype focused
phenoCentilePlots <- c()

phenoCentilePlots[[1]] <- ggplot(data=violinDf, aes(regions, centilesList)) +
  geom_violin(color="gray", fill="gray", alpha=0.5) +
  geom_jitter(height = 0, width=0.15, aes(color=scanYearGroup), alpha=0.65) +
  scale_color_manual(values = cbbPalette, name = "Year of Scan") +
  labs(title="Distributions of Centiles (FS)") + 
  xlab("Tissue Type") +
  ylab("Centile") + 
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

for (r in c(1:length(phenos))){
  phenoDf <- violinDf[regions==phenos[[r]], ]
  phenoCentilePlots[[r+1]] <- ggplot(data=phenoDf, aes(scanYearGroup, centilesList)) +
    geom_violin(color="gray", fill="gray", alpha=0.35) +
    geom_jitter(height = 0, width=0.15, aes(color=scanYearGroup), alpha=0.65) +
    scale_color_manual(values = cbbPalette, name = "Year of Scan") +
    labs(title=paste0("Centile (Phenotype = ", phenos[[r]],")")) + 
    xlab("Scanner ID") +
    ylab("Centile") + 
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 18))
}
patch <- wrap_plots(phenoCentilePlots, ncol=2, guides="collect")
png(file="/Users/youngjm/Data/slip/figures/2022-11-08_clip_predicted_centiles_scanners.png",
    width=1200, height=1000)
print(patch)
dev.off()

write.csv(ssDf, "/Users/youngjm/Data/slip/fs6_stats/07_fully_filtered_postcombat_clip_ss.csv")

