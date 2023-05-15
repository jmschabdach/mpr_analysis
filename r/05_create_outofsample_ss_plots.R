gc()
setwd("/Users/youngjm/Projects/mpr_analysis/r/")
source("lib_mpr_analysis.r")

#-------------------------------------------------------------------------------
# Filenames to update
#-------------------------------------------------------------------------------

# Load the dataframe containing subject demographics, imaging phenotypes, etc.
baseDir <- '/Users/youngjm/Data/slip/fs6_stats/'
inFn <- paste0(baseDir, '2023-04-24_supplemental_synthseg_phenotypes.csv')
demoFn <- paste0(baseDir, '2023-04-24_supplemental_demo.csv')
fnOut <- paste0(baseDir, '2023-04-24_supplemental_synthseg_phenotypes.csv')
scannerFnOut <- paste0(baseDir, 'ss_scanner_id_lookup_table.csv')
fnLbcc <- '/Users/youngjm/Data/lifespan_growth_charts/slip-agemedian_centiles.csv'
fnSlipSs <- paste0(baseDir, "07_fully_filtered_postcombat_clip_ss.csv")

#-------------------------------------------------------------------------------
# Load and clean data
#-------------------------------------------------------------------------------


demoDf <- read.csv(demoFn)
outOfSampleDf <- read.csv(inFn)

demoDf$subj_id <- paste0("sub-", demoDf$pat_id)
# Drop columns we don't need any more (nf columns)
toDrop <- c("X", "pat_id", "proc_ord_id", "age_in_days", "arcus_id") 
demoDf <- demoDf[ , -which(names(demoDf) %in% toDrop)]
# Map Female to F and Male to M
demoDf <- demoDf %>%
  mutate(sex = case_when(
    grepl('Female', sex, fixed=TRUE) ~ 'F',
    grepl('Male', sex, fixed=TRUE) ~ 'M'
  ))

outOfSampleDf <- merge(outOfSampleDf, demoDf, by="subj_id")

if ('age_in_days' %in% colnames(outOfSampleDf)){
  names(outOfSampleDf)[names(outOfSampleDf) == 'age_in_days'] <- 'age_at_scan_days'
}

if ('subj_id' %in% colnames(outOfSampleDf)){
  names(outOfSampleDf)[names(outOfSampleDf) == 'subj_id'] <- 'patient_id'
}

analysisDf <- outOfSampleDf
analysisDf$patient_id <- droplevels(as.factor(analysisDf$patient_id))
print(length(levels(analysisDf$patient_id)))

# Make an age in years column from the age in days column
analysisDf$age_in_years <- (analysisDf$age_at_scan_days)/365.25

# Some of the columns in this df should be factors
toFactor <- c('sex', 'fs_version')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

# We only one scan per subject
# Sort the dataframe by patient_id and scanner_id
analysisDf <- analysisDf[ with(analysisDf, order(analysisDf$patient_id, analysisDf$scan_id)), ]
# Drop all but the first occurrence of each patient_id
analysisDf <- analysisDf[!duplicated(analysisDf$patient_id), ]
# Convert patient_id to a factor - idk why I did this?
analysisDf$patient_id <- droplevels(as.factor(analysisDf$patient_id))

if ('CortexVol' %in% colnames(analysisDf)){
  analysisDf$GMV <- analysisDf$CortexVol
  analysisDf$WMV <- analysisDf$CerebralWhiteMatterVol
  analysisDf$sGMV <- analysisDf$SubCortGrayVol
  analysisDf$SA <- analysisDf$CorticalSurfaceArea
  analysisDf$CT <- analysisDf$MeanCorticalThickness
} else if ('Cortex' %in% colnames(analysisDf)){
  analysisDf$GMV <- analysisDf$Cortex
}

analysisDf$TCV <- analysisDf$GMV + analysisDf$WMV


if ('BrainSegVol' %in% colnames(analysisDf)){
  analysisDf$CSF <- analysisDf$BrainSegVol - analysisDf$BrainSegVolNotVent
} else if ('Ventricles' %in% colnames(analysisDf)) {
  analysisDf$CSF <- analysisDf$Ventricles
}

# Drop any scans with NAs
analysisDf <- analysisDf[complete.cases(analysisDf), ]

# Drop any scans with bad phenotype values
analysisDf[, "CSF"][analysisDf[, "CSF"] <= 0] <- NA


## -----------------------------------------------------------------------------

# Set up variables needed for the analysis
phenos <- c('GMV', 'WMV', 'sGMV', 'CSF', 'TCV') 
cbbPalette <- c("#E69F00", "#000000")

lbccDf <- read.csv(fnLbcc)
slipSsDf <- read.csv(fnSlipSs) 

# SLIP age lc
lbccDf$feature <- as.factor(lbccDf$feature)
names(lbccDf)[names(lbccDf) == "age"] <- "logAge"
lbccDf$age <- (10^lbccDf$logAge) # post conception age in days

# Convert scanner id to a factor
slipSsDf$scanner_id <- as.factor(slipSsDf$scanner_id)
slipSsDf$sex <- as.factor(slipSsDf$sex)
slipSsDf$proc_ord_year <- as.factor(slipSsDf$proc_ord_year)
slipSsDf$top_scan_reason_factors <- as.factor(slipSsDf$top_scan_reason_factors)
# # Replace any 0 valued phenotypes with NA
# slipSsDf[, "WMV"][slipSsDf[, "WMV"] <= 0] <- NA
# slipSsDf[, "CSF"][slipSsDf[, "CSF"] <= 0] <- NA

slipSsDf <- na.omit(slipSsDf)
analysisDf <- na.omit(analysisDf)

# Determine the limited age range for predicting the centiles with the GAMLSS
analysisDf$logAge <- log(analysisDf$age_at_scan_days+280, 10)
minAgeDaysPostCon <- min(min(analysisDf$logAge), min(slipSsDf$logAge))
maxAgeDaysPostCon <- max(max(analysisDf$logAge), max(slipSsDf$logAge))

ageLimited <- unique(lbccDf[lbccDf$logAge > minAgeDaysPostCon & lbccDf$logAge < maxAgeDaysPostCon, "logAge"])
lbccDf <- lbccDf[lbccDf$logAge > minAgeDaysPostCon & lbccDf$logAge < maxAgeDaysPostCon, ]

# Set up a list of tick marks to use on log(post-conception age) x-axes
tickMarks <- c()
for (year in c(0, 1, 2, 5, 10, 20)){ # years
  tickMarks <- append(tickMarks, log(year*365.25 + 280, base=10))
}
tickLabels <- c("Birth", "1", "2", "5", "10", "20")

slipSsDf$stage <- as.factor("Primary")
analysisDf$stage <- as.factor("Out of Sample")

#-------------------------------------------------------------------------------
# Part 2: Growth charts for global phenotypes
#-------------------------------------------------------------------------------

# Set up lists for storing variables in the following for loops
centilesOrig <- c()
centilesOos <- c()

for (p in phenos) {
  print(p)
  smallLc <- lbccDf[lbccDf$feature==p, ]
  plt <- c()
  
  # 1. Build the models for the primary and out of sample data
  formulaSs <- as.formula(paste0(p, "~fp(logAge, npoly=3) + sex - 1"))
  gamModelSs <-gamlss(formula = formulaSs,
                      sigma.formula = formulaSs,
                      nu.formula = as.formula(paste0(p, "~1")),
                      family = GG,
                      data = na.omit(slipSsDf),
                      # data = rbind(na.omit(clipDf), placeholderDf),
                      control = gamlss.control(n.cyc = 200),  # lifespan
                      trace = F)
  predictedMedianOrigSs <- predictCentilesForAgeRange(gamModelSs, ageLimited)

  fanCentiles <- c()
  desiredCentiles <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
  for (i in c(1:length(desiredCentiles))){
    fanCentiles[[i]] <- predictCentilesForAgeRange(gamModelSs, ageLimited, cent=desiredCentiles[[i]])
  }
  
  centilesOrig <- calculatePhenotypeCentile(gamModelSs, slipSsDf[[p]], slipSsDf$logAge, slipSsDf$sex)
  centilesOos <- calculatePhenotypeCentile(gamModelSs, analysisDf[[p]], analysisDf$logAge, analysisDf$sex)
  
  # Get the y min and y max values
  ymin <- 0.85*min(min(analysisDf[,p]), min(slipSsDf[,p]))
  ymax <- 1.05*max(max(analysisDf[,p]), max(slipSsDf[,p]))
  xmax <- max(max(analysisDf$logAge), max(slipSsDf$logAge))
  
  # 5. Plot CLIP vs Lifespan
  if (p == "CSF") {
    ageLimited[which(fanCentiles[[1]] < 50000)]
    
    ymin <- 0.85*min(min(analysisDf[,p]), min(slipSsDf[,p]))
    ymax <- max(1.05*max(analysisDf[,p]), 1.05*max(slipSsDf[,p]), 55000)
    
    plt[[1]] <- ggplot(x="linear") +
      geom_point(data=analysisDf, aes(x=logAge, y=analysisDf[,p], color=stage), alpha=0.3) +
      geom_point(data=slipSsDf, aes(x=logAge, y=slipSsDf[,p], color=stage), alpha=0.3) +
      geom_line(aes(x=ageLimited[which(fanCentiles[[1]] < 50000)], y=fanCentiles[[1]][which(fanCentiles[[1]] < 50000)], linetype=" 5th")) +
      geom_line(aes(x=ageLimited[which(fanCentiles[[2]] < 50000)], y=fanCentiles[[2]][which(fanCentiles[[2]] < 50000)], linetype="10th")) +
      geom_line(aes(x=ageLimited[which(fanCentiles[[3]] < 50000)], y=fanCentiles[[3]][which(fanCentiles[[3]] < 50000)], linetype="25th")) +
      geom_line(aes(x=ageLimited[which(fanCentiles[[4]] < 50000)], y=fanCentiles[[4]][which(fanCentiles[[4]] < 50000)], linetype="50th")) +
      geom_line(aes(x=ageLimited[which(fanCentiles[[5]] < 50000)], y=fanCentiles[[5]][which(fanCentiles[[5]] < 50000)], linetype="75th")) +
      geom_line(aes(x=ageLimited[which(fanCentiles[[6]] < 50000)], y=fanCentiles[[6]][which(fanCentiles[[6]] < 50000)], linetype="90th")) +
      geom_line(aes(x=ageLimited[which(fanCentiles[[7]] < 50000)], y=fanCentiles[[7]][which(fanCentiles[[7]] < 50000)], linetype="95th")) +
      scale_color_manual(values = cbbPalette, name = "Stage of Data") +
      # labs(subtitle = paste0(p," SynthSeg+ (r(orig, new)=", format(rSsOrigNew[[p]],digits=3),")")) +
      xlab("Age at scan (log(years))") +
      scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                         limits=c(tickMarks[[1]], xmax)) +
      scale_discrete_manual("linetype", values=c(" 5th"='dotted', 
                                                 "10th"='dotdash', 
                                                 "25th"='dashed', 
                                                 "50th"='solid', 
                                                 "75th"='dashed', 
                                                 "90th"='dotdash', 
                                                 "95th"='dotted')) +
      labs(linetype="Growth Centile Curves") +
      ylab(paste0(p, " Volume")) + 
      ylim(ymin, ymax) +
      theme(axis.line = element_line(colour = "black"),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 14))
  } else {
    fanCat <- c(" 5th", "10th", "25th", "50th", "75th", "90th", "95th")
    plt[[1]] <- ggplot(x="linear") +
      geom_point(data=analysisDf, aes(x=logAge, y=analysisDf[,p], color=stage), alpha=0.5) +
      geom_point(data=slipSsDf, aes(x=logAge, y=slipSsDf[,p], color=stage), alpha=0.3) +
      geom_line(aes(x=ageLimited, y=fanCentiles[[1]], linetype=" 5th")) +
      geom_line(aes(x=ageLimited, y=fanCentiles[[2]], linetype="10th")) +
      geom_line(aes(x=ageLimited, y=fanCentiles[[3]], linetype="25th")) +
      geom_line(aes(x=ageLimited, y=fanCentiles[[4]], linetype="50th")) +
      geom_line(aes(x=ageLimited, y=fanCentiles[[5]], linetype="75th")) +
      geom_line(aes(x=ageLimited, y=fanCentiles[[6]], linetype="90th")) +
      geom_line(aes(x=ageLimited, y=fanCentiles[[7]], linetype="95th")) +
      scale_color_manual(values = cbbPalette, name = "Stage of Data") +
      xlab("Age at scan (log(years))") +
      scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                         limits=c(tickMarks[[1]], xmax)) +
      scale_discrete_manual("linetype", values=c(" 5th"='dotted', 
                                     "10th"='dotdash', 
                                     "25th"='dashed', 
                                     "50th"='solid', 
                                     "75th"='dashed', 
                                     "90th"='dotdash', 
                                     "95th"='dotted')) +
      labs(linetype="Growth Centile Curves") +
      ylab(paste0(p, " Volume")) + 
      ylim(ymin, ymax) +
      theme(axis.line = element_line(colour = "black"),
            # panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 14))
  }
  
  stage <- c(rep('Primary', length(centilesOrig)), rep('Out of Sample', length(centilesOos)))
  centilesList <- c(centilesOrig, centilesOos)
  violinDf <- data.frame(stage, centilesList)
  
  stat.test <- t_test(violinDf, as.formula("centilesList ~ stage"))
  testResults <- ks.test(analysisDf[,p], slipSsDf[,p])
  stat.test[, p] <- as.numeric(testResults$p.value)
  testResults <- ks.test(analysisDf[,p], slipSsDf[,p])
  if (testResults$p.value < (0.05/6)) {
    print("Statistically significant DIFFERENCE")
  } else {
    print("No difference in distributions")
    stat.test$p.adj.signif <- c("ns")
  }
  
  # Generate violin plots for primary and out of sample data
  plt[[2]] <- ggplot(data=violinDf, aes(stage, centilesList)) +
    geom_violin(color="gray", fill="gray", alpha=0.5) +
    geom_jitter(height = 0, width=0.15, aes(color=stage), alpha=0.65) +
    stat_pvalue_manual(stat.test, y.position= 1.05, label = "p.adj.signif") +
    scale_color_manual(values = cbbPalette, labels = c("Primary", "Out of Sample"), name = "Dataset") +
    xlab(p) +
    ylab("Centile") + 
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          text = element_text(size = 14))
  
  # Build the plot
  design <- "AAB"
  patch <- wrap_plots(plt, nrow = 1, design = design, guides='collect')
  png(file=paste0("/Users/youngjm/Data/slip/figures/2023_newdata_",p,".png"),
      width=900, height=400)
  print(patch) 
  dev.off()
}

##------------------------------------------------------------------------------
# Get info for demographics table

# Pull the stats for the table
print(table(analysisDf$sex))
analysisDf$age_group <- case_when(
  analysisDf$age_in_years < 2.0 ~ "0-2",
  analysisDf$age_in_years < 5.0 ~ "02-5",
  analysisDf$age_in_years < 10.0 ~ "05-10",
  analysisDf$age_in_years < 13.0 ~ "10-13",
  analysisDf$age_in_years < 18.0 ~ "13-18",
  TRUE ~ "18+"
)
analysisDf$age_group <- as.factor(analysisDf$age_group)

analysisDf$scan_year_group <- case_when(
  analysisDf$proc_ord_year %in% c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020") ~ "2011+",
  TRUE ~ as.character(analysisDf$proc_ord_year)
)
analysisDf$scan_year_group <- as.factor(analysisDf$scan_year_group)

print(table(analysisDf$age_group))
print(table(analysisDf[analysisDf$sex == "M", "age_group"]))
print(table(analysisDf[analysisDf$sex == "F", "age_group"]))

print(table(analysisDf$scan_year_group))
print(table(analysisDf[analysisDf$sex == "M", "scan_year_group"]))
print(table(analysisDf[analysisDf$sex == "F", "scan_year_group"]))

analysisDf$race <- as.factor(analysisDf$race)
print(table(analysisDf$race))
print(table(analysisDf[analysisDf$sex == "M", "race"]))
print(table(analysisDf[analysisDf$sex == "F", "race"]))

