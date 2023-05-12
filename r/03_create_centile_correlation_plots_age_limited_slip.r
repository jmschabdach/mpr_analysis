gc()
setwd("/Users/youngjm/Projects/mpr_analysis/r/")
source("lib_mpr_analysis.r")

##
# Note:
# Please change the paths in the following few lines to match your directory 
# structure. 

#-------------------------------------------------------------------------------
# Load the data from multiple files and make sure all of the data types match
#-------------------------------------------------------------------------------
# CSV input file names
fn <- '/Users/youngjm/Data/lifespan_growth_charts/slip-agemedian_centiles.csv'
lbccDf <- read.csv(fn)
fn <- '/Users/youngjm/Data/slip/fs6_stats/06_combatted_fs_plus_metadata.csv'
slipFsDf<- read.csv(fn)
fn <- '/Users/youngjm/Data/slip/fs6_stats/06_combatted_ss_plus_metadata.csv'
slipSsDf <- read.csv(fn) 
preCombatFn <- '/Users/youngjm/Data/slip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
preCombatDf <- read.csv(preCombatFn)

# Figure output file names
figFnPredictedCentiles <- "/Users/youngjm/Data/slip/figures/2022-11-08_clip_predicted_centiles.png"
figFnPredictedCentilesByScanner <- "/Users/youngjm/Data/slip/figures/2022-11-08_clip_predicted_centiles_scanners.png"
figFnPhenoCorrelationBase <- "/Users/youngjm/Data/slip/figures/2022-10-19_clip_agelimited_lifespan_correlation_"

# CSV output file names
fnFsCentiles <- "/Users/youngjm/Data/slip/fs6_stats/fs_gamlss_centiles.csv"
fnOutFilteredFsSlip <- "/Users/youngjm/Data/slip/fs6_stats/07_fully_filtered_postcombat_clip_fs.csv"
fnOutFilteredSsSlip <- "/Users/youngjm/Data/slip/fs6_stats/07_fully_filtered_postcombat_clip_ss.csv"


# Prep synthseg data
# Add a column to synthseg for post-conception age
slipSsDf$age <- slipSsDf$age_at_scan_days+280
# Add log of post-conception age
slipSsDf$logAge <- log(slipSsDf$age, 10)
# Convert scanner id to a factor
slipSsDf$scanner_id <- gsub("Scanner ","",as.character(slipSsDf$scanner_id))
slipSsDf$scanner_id <- as.factor(slipSsDf$scanner_id)
slipSsDf$sex <- as.factor(slipSsDf$sex)
slipSsDf$proc_ord_year <- as.factor(slipSsDf$proc_ord_year)
slipSsDf$top_scan_reason_factors <- as.factor(slipSsDf$top_scan_reason_factors)
# Rename columns to match the slipFsDfcolumns
names(slipSsDf)[names(slipSsDf) == "Cortex"] <- 'GMV'
names(slipSsDf)[names(slipSsDf) == "Ventricles"] <- 'CSF'
names(slipSsDf)[names(slipSsDf) == "subj_id"] <- 'patient_id'
# Replace any 0 valued phenotypes with NA
slipSsDf[, "WMV"][slipSsDf[, "WMV"] <= 0] <- NA
slipSsDf[, "CSF"][slipSsDf[, "CSF"] <= 0] <- NA

# Prep CLIP FS data
# Add a column to clip for log(post conception age)
slipFsDf$logAge <- log(slipFsDf$age_at_scan_days+280, 10)
# Make sex and scanner_id factors
slipFsDf$sex <- as.factor(slipFsDf$sex)
slipFsDf$scanner_id <-gsub("Scanner ","",as.character(slipFsDf$scanner_id))
slipFsDf$scanner_id <- as.factor(slipFsDf$scanner_id)
slipFsDf$proc_ord_year <- as.factor(slipFsDf$proc_ord_year)
slipFsDf$top_scan_reason_factors <- as.factor(slipFsDf$top_scan_reason_factors)
# Drop rows with missing values
slipFsDf<- drop_na(slipFsDf)

# Prep the precombat FS data
preCombatDf$sex <- as.factor(preCombatDf$sex)
preCombatDf$logAge <- log((preCombatDf$age_at_scan_days + 280), base=10)
preCombatDf$scanner_id <- gsub("Scanner ","",as.character(preCombatDf$scanner_id))
preCombatDf$scanner_id <- as.factor(preCombatDf$scanner_id)
preCombatDf$top_scan_reason_factors <- as.factor(preCombatDf$top_scan_reason_factors)
preCombatDf$CSF <- abs(preCombatDf$BrainSegVol - preCombatDf$BrainSegVolNotVent)
preCombatDf <- drop_na(preCombatDf)

# Do an inner join on the pre and post combat FS dataframes
slipFsDf<- inner_join(slipFsDf, preCombatDf,
           by=c("patient_id", "scanner_id", "age_at_scan_days", "logAge", "top_scan_reason_factors"),
           suffix=c("", ".pre"))

# Filter out any subjects not in both FS and SS dataframes
slipSsDf <- drop_na(slipSsDf)
slipSsDf <- slipSsDf[(slipSsDf$patient_id %in% slipFsDf$patient_id), ]
slipFsDf<- slipFsDf[(slipFsDf$patient_id %in% slipSsDf$patient_id), ]

# Order the data by age at scan
setorder(slipFsDf, age_at_scan_days)  
setorder(slipSsDf, age_at_scan_days)

# Set up variables needed for the analysis
phenos <- c('GMV', 'WMV', 'sGMV', 'CSF', 'TCV')  
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",  
                "#D55E00", "#CC79A7", "#0072B2", "#F0E442")

# Determine the limited age range for predicting the centiles with the GAMLSS
lbccDf$feature <- as.factor(lbccDf$feature)
names(lbccDf)[names(lbccDf) == "age"] <- "logAge"
lbccDf$age <- (10^lbccDf$logAge) # post conception age in days

minAgeDaysPostCon <- min(slipFsDf$logAge)
maxAgeDaysPostCon <- max(slipFsDf$logAge)
ageLimited <- unique(lbccDf[lbccDf$logAge > minAgeDaysPostCon & lbccDf$logAge < maxAgeDaysPostCon, "logAge"])
lbccDf <- lbccDf[lbccDf$logAge > minAgeDaysPostCon & lbccDf$logAge < maxAgeDaysPostCon, ]

#-------------------------------------------------------------------------------
# Part 2: Growth charts for global phenotypes
#-------------------------------------------------------------------------------

# Set up lists for storing variables in the following for loops
rFsLbcc <- c()
rSsLbcc <- c()
rFsSs <- c()
rFsSsPhenos <- c()

ageAtPeakFs <- c()
ageAtPeakSs <- c()
ageAtPeakLbcc <- c()

centilesFs <- c()
centilesFsPreCombat <- c()
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
  fnOut <- paste0(figFnPhenoCorrelationBase, p, '.png')

  ## Generate GAMLSS models
  # FreeSurfer model: uses SurfaceHoles (euler number)
  formulaFs <- as.formula(paste0(p, "~fp(logAge, npoly=3) + SurfaceHoles + sex - 1"))
  gamModelFs <-gamlss(formula = formulaFs,
                    sigma.formula = formulaFs,
                    nu.formula = as.formula(paste0(p, "~1")),
                    family = GG,
                    data = na.omit(slipFsDf),
                    control = gamlss.control(n.cyc = 200),  # lifespan
                    trace = F)
  
  # SynthSeg model: does not produce SurfaceHoles, so they are removed
  formulaSs <- as.formula(paste0(p, "~fp(logAge, npoly=3) + sex - 1"))
  gamModelSs <-gamlss(formula = formulaSs,
                      sigma.formula = formulaSs,
                      nu.formula = as.formula(paste0(p, "~1")),
                      family = GG,
                      data = na.omit(slipSsDf),
                      control = gamlss.control(n.cyc = 200),  # lifespan
                      trace = F)
  print("finished training the models")

  ## Predict phenotype values for set age range
  predictedMedianFs <- predictCentilesForAgeRange(gamModelFs, ageLimited, median(slipFsDf$SurfaceHoles))
  predictedMedianSs <- predictCentilesForAgeRange(gamModelSs, ageLimited)
  
  
  # 3. Get a subset of the Lifespan dataframe specific to the current phenotype
  lbccPheno <- lbccDf[lbccDf$feature==p, ]
  
  ## Calculate various statistics
  # Correlation 
  rFsLbcc[[p]] <- cor(predictedMedianFs, lbccPheno$value)
  rSsLbcc[[p]] <- cor(predictedMedianSs, lbccPheno$value)
  rFsSsPhenos[[p]] <- cor(slipFsDf[,p], slipSsDf[,p])
  
  # Age at peak phenotype value
  ageAtPeakFs[[p]] <- 10^(sort(ageLimited)[which.max(predictedMedianFs)])
  ageAtPeakSs[[p]] <- 10^(sort(ageLimited)[which.max(predictedMedianSs)])
  ageAtPeakLbcc[[p]] <- lbccPheno[lbccPheno$value == max(lbccPheno$value), 'age']
  
  ## Grab data to use later in violin plots
  centilesFs[[p]] <- calculatePhenotypeCentile(gamModelFs, slipFsDf[[p]], slipFsDf$logAge, slipFsDf$sex, slipFsDf$SurfaceHoles)
  centilesFsPreCombat[[p]] <- calculatePhenotypeCentile(gamModelFs, slipFsDf[[paste0(p, ".pre")]], slipFsDf$logAge, slipFsDf$sex, slipFsDf$SurfaceHoles)
  regions <- append(regions, rep(p, length(centilesFs[[p]])))
  idxes <- append(idxes, c(1:length(centilesFs[[p]])))
  reasons <- append(reasons, slipFsDf$top_scan_reason_factors)
  yearOfScan <- append(yearOfScan, slipFsDf$proc_ord_year)
  
  # Get the y min and y max values
  ymin <- min(min(slipFsDf[,p]), min(slipSsDf[,p]))
  ymax <- max(max(slipFsDf[,p]), max(slipSsDf[,p]))
  
  # 5. Plot CLIP vs Lifespan
  plots[[1]] <- ggplot() +
    geom_point(data=slipFsDf, aes(x=logAge, y=slipFsDf[,p], color=top_scan_reason_factors), alpha=0.3) +
    geom_line(aes(x=ageLimited, y=predictedMedianFs, linetype="Predicted for SLIP")) +
    geom_line(data=lbccPheno, aes(x=logAge, y=value, linetype="SLIP-age LBCC")) +
    scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
    scale_linetype_manual(values = c('solid', 'dashed'), name="50th Centile")+
    scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                       limits=c(tickMarks[[1]], max(slipFsDf$logAge))) +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(subtitle = paste0("FreeSurfer v. 6.0.0 (r=", format(rFsLbcc[[p]], digits=3), ")")) +
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
    geom_point(data=slipSsDf, aes(x=logAge, y=slipSsDf[,p], color=top_scan_reason_factors), alpha=0.3) +
    geom_line(aes(x=ageLimited, y=predictedMedianSs, linetype="Predicted for SLIP")) +
    geom_line(data=lbccPheno, aes(x=logAge, y=value, linetype="SLIP-age LBCC")) +
    scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
    scale_linetype_manual(values = c('solid', 'dashed'), name="50th Centile")+
    labs(subtitle = paste0("SynthSeg+ (r=", format(rSsLbcc[[p]], digits=3), ")")) +
    xlab("Age at scan (log(years))") +
    scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                       limits=c(tickMarks[[1]], max(slipFsDf$logAge))) +
    ylab("Phenotype Value") + 
    ylim(ymin, ymax) +
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 18))
  
  plots[[3]] <- ggplot() +
    geom_point(aes(x=slipFsDf[,p], y=slipSsDf[,p], color=slipFsDf$top_scan_reason_factors), alpha=0.3) +
    scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
    labs(subtitle = paste0("FreeSurfer vs SynthSeg (r=", format(rFsSsPhenos[[p]], digits=3), ")")) +
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


#-------------------------------------------------------------------------------
# Generate a plot with the centile fan and the centile violin plots
#-------------------------------------------------------------------------------

# Generate GAMLSS models for FreeSurfer GMV phenotype
formulaGmv <- as.formula("GMV~fp(logAge, npoly=3) + SurfaceHoles + sex - 1")
gmvFsModel <-gamlss(formula = formulaGmv,
                  sigma.formula = formulaGmv,
                  nu.formula = as.formula("GMV~1"),
                  family = GG,
                  data = na.omit(slipFsDf),
                  control = gamlss.control(n.cyc = 200),  # lifespan
                  trace = F)
print("finished training the models")

# Prep list for storing calculated centiles
fanCentiles <- c()

# For each desired centile, calculate the centile values for the FreeSurfer GMV model
desiredCentiles <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
for (i in c(1:length(desiredCentiles))){
  fanCentiles[[i]] <- predictCentilesForAgeRange(gmvFsModel, ageLimited, 
                                                 median(slipFsDf$SurfaceHoles), 
                                                 desiredCentiles[[i]])
}

## Make an age vs. phenotype plot of the various centile lines
sampleCentileFan <- ggplot() +
  geom_point(aes(x=slipFsDf$logAge, slipFsDf$GMV), alpha=0.5) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[1]]), alpha=0.4) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[2]]), alpha=0.6) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[3]]), alpha=0.8) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[4]])) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[5]]), alpha=0.8) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[6]]), alpha=0.6) +
  geom_line(aes(x=ageLimited, y=fanCentiles[[7]]), alpha=0.4) +
  # scale_color_manual(values = cbbPalette, name = "Phenotype Value") +
  scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(slipFsDf$logAge))) +
  labs(title="Sample Centile Growth Chart for GMV") + 
  xlab("Age at scan (log(years))") +
  ylab("GMV Centile") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

# Make age vs. phenotype scatterplot
rawPhenoPlot <- ggplot(x="linear") +
  geom_point(data=slipFsDf, aes(x=logAge, y=GMV, color=top_scan_reason_factors), alpha=0.3) +
  scale_color_manual(values = cbbPalette, labels = reasonLabels, name = "Reason for Scan") +
  labs(title="Phenotype Distribution for GMV") + 
  xlab("Age at scan (log(years))") +
  ylab("GMV Value") + 
  scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(slipFsDf$logAge))) +
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 18))

## Make violin plots showing the distribution of phenotype centiles across reason for scan
centilesFsList <- c(centilesFs$GMV, centilesFs$WMV, centilesFs$sGMV, centilesFs$CSF, centilesFs$TCV)
vScanId <- c(rep(slipFsDf$scan_id, 5))
violinDf <- data.frame(idxes, regions, centilesFsList, reasons, yearOfScan, vScanId)
write.csv(violinDf, fnFsCentiles, row.names = FALSE)

violin <- ggplot(data=violinDf, aes(regions, centilesFsList)) +
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


## Make violin plots showing the distribution of phenotype centiles across reason for scan
phenoCentilePlots <- c()
phenoCentilePlots[[1]] <- violin
for (p in c(1:length(phenos))){
  phenoDf <- violinDf[regions==phenos[[p]], ]
  phenoCentilePlots[[p+1]] <- ggplot(data=phenoDf, aes(reasons, centilesFsList)) +
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

# Combine plots from the previous 
topPanel <- wrap_plots(rawPhenoPlot + sampleCentileFan)
bottomPanel <- wrap_plots(phenoCentilePlots, ncol=2, guides="collect")
layout <-"
A
B
B
B
"
patch <- wrap_plots(topPanel + bottomPanel + plot_layout(design=layout), guides="auto")
png(file=figFnPredictedCentiles,
    width=1400, height=1400)
print(patch)
dev.off()

#-------------------------------------------------------------------------------
# Create violin plots showing the distribution of phenotype centiles among scanners
#-------------------------------------------------------------------------------
scanYearGroup <- case_when(
  yearOfScan %in% c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020") ~ "2011+",
  TRUE ~ as.character(yearOfScan)
)
violinDf <- data.frame(idxes, regions, centilesFsList, reasons, scanYearGroup)

phenoCentilePlots <- c()
phenoCentilePlots[[1]] <- ggplot(data=violinDf, aes(regions, centilesFsList)) +
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
  phenoCentilePlots[[r+1]] <- ggplot(data=phenoDf, aes(scanYearGroup, centilesFsList)) +
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
png(file=figFnPredictedCentilesByScanner,
    width=1200, height=1000)
print(patch)
dev.off()


#-------------------------------------------------------------------------------
# Statistical test: are the distributions of centiles before and after ComBat comparable?
#-------------------------------------------------------------------------------

# Test the distributions of phenotype values before and after ComBat
plots <- c()
pvals <- c()
for (i in seq(1:length(phenos))){
  print(i)
  p <- phenos[[i]]
  # Set up the data frame
  tmp <- data.frame(PreCombat= centilesFsPreCombat[[i]],
                    PostCombat = centilesFs[[i]],
                    scannerId = slipFsDf$scanner_id)
  reps <- length(tmp$PreCombat)
  
  # Transform the data
  tmp <- data.frame(Time = c(rep("pre-ComBat", reps), rep("post-ComBat", reps)),
                    Phenotype=c(tmp$PreCombat, tmp$PostCombat),
                    Scanner = c(tmp$scannerId, tmp$scannerId))
  
  # Reorder Time
  tmp$Time <- factor(tmp$Time, levels = c("pre-ComBat", "post-ComBat"))
  tmp$Scanner <- factor(tmp$Scanner)
  
  # Make plots for pre- and post-Combat FreeSurfer phenotype centiles
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
    # Make the plots
    plots[[paste0(p, "_", time)]] <- ggplot(tmp[tmp$Time == time, ]) +
      aes(x=Scanner, y=Phenotype) +
      geom_boxplot(aes(fill=Scanner), varwidth = TRUE, alpha=0.75) +
      geom_jitter(height = 0, width=0.15, alpha=0.4) +
      scale_fill_viridis_d(name = "Scanner ID") +
      labs(title=paste0(p, " ", time),
           subtitle=paste0("(F(",paste(dim(slipFsDf)[1]-1), ", ", paste(dof),") = ",
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

# Save the plots
imgOut <- "~/Data/slip/figures/2022-11-16_pre_post_combat_combo_all.png"
patch <- wrap_plots(plots, ncol=2, byrow=TRUE, guides="collect")
png(file=imgOut,
    width=1000, height=1400)
print(patch)
dev.off()

# Statistical tests: for each phenotype, compare the centile values across reason for scan (ANOVA)
# centile ~ reason for scan | phenotype
violinDf$reasons <- as.factor(violinDf$reasons)
violinDf$yearOfScan <- as.factor(yearOfScan)
for (p in phenos){
  aovOut <- aov(centilesFsList ~ yearOfScan, data=na.omit(violinDf[regions == p,]))
  # lm , yearOfScan = numeric
  print(summary(aovOut))
}

#-------------------------------------------------------------------------------
# Print info for tables
#-------------------------------------------------------------------------------

convertAgeToYears(ageAtPeakFs)
convertAgeToYears(ageAtPeakSs)
convertAgeToYears(ageAtPeakLbcc)

# Pull the stats for the table
print(table(slipFsDf$sex))
slipFsDf$age_group <- case_when(
  slipFsDf$age_in_years < 2.0 ~ "0-2",
  slipFsDf$age_in_years < 5.0 ~ "02-5",
  slipFsDf$age_in_years < 10.0 ~ "05-10",
  slipFsDf$age_in_years < 13.0 ~ "10-13",
  slipFsDf$age_in_years < 18.0 ~ "13-18",
  TRUE ~ "18+"
)
slipFsDf$age_group <- as.factor(slipFsDf$age_group)

slipFsDf$scan_year_group <- case_when(
  slipFsDf$proc_ord_year %in% c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020") ~ "2011+",
  TRUE ~ as.character(slipFsDf$proc_ord_year)
)
slipFsDf$scan_year_group <- as.factor(slipFsDf$scan_year_group)

print(table(slipFsDf$age_group))
print(table(slipFsDf[slipFsDf$sex == "M", "age_group"]))
print(table(slipFsDf[slipFsDf$sex == "F", "age_group"]))

print(table(slipFsDf$top_scan_reason_factors))
print(table(slipFsDf[slipFsDf$sex == "M", "top_scan_reason_factors"]))
print(table(slipFsDf[slipFsDf$sex == "F", "top_scan_reason_factors"]))

print(table(slipFsDf$scan_year_group))
print(table(slipFsDf[slipFsDf$sex == "M", "scan_year_group"]))
print(table(slipFsDf[slipFsDf$sex == "F", "scan_year_group"]))

write.csv(slipFsDf, fnOutFilteredFsSlip)
write.csv(slipSsDf, fnOutFilteredSsSlip)
