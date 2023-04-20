library(ggplot2)
library(gamlss)
library(tidyr)
library(dplyr)

fnBase <- "/Users/youngjm/Data/lifespan_growth_charts/"
ogFn <- "/Users/youngjm/Data/lifespan_growth_charts/lifespan_subset.csv"
ogDf <- read.csv(ogFn)

print(dim(ogDf))
regionCols <- c("GMV", "WMV", "sGMV", "CSF")
metaCols <- c("scan_id", "sex", "study", "fs_version", "site", "session",
              "age_days") #, "new_euler")
keepCols <- append(metaCols, regionCols)
names(ogDf)[names(ogDf) == "participant"] = "scan_id"
names(ogDf)[names(ogDf) == "Ventricles"] = "CSF"

# Keep only the columns we care about
analysisDf <- subset(ogDf, select=keepCols)
# Replace new_euler NA with median
# analysisDf$new_euler[is.na(analysisDf$new_euler)]<-median(analysisDf$new_euler,na.rm=TRUE)
analysisDf <- drop_na(analysisDf)
analysisDf$site <- paste0("prefix_", analysisDf$site)


# Filter for only FS6_T1 and FSInfant processing
analysisDf <- analysisDf[(analysisDf$fs_version == "FS6_T1" | analysisDf$fs_version == "FSInfant"),]

# Filter by age
ageCutoff <- (365.25*22) + 280
analysisDf <- analysisDf[analysisDf$age_days <= ageCutoff, ]

# Convert sex, fs_version, site all to factors
analysisDf$sex <- as.factor(analysisDf$sex)
analysisDf$site <- as.factor(analysisDf$site)
analysisDf$fs_version <- as.factor(analysisDf$fs_version)

# Add site count
tmp <- table(analysisDf$site)
analysisDf <- mutate(analysisDf, site_popularity=tmp[site])
# Drop sites occurring fewer than X times
analysisDf <- analysisDf[analysisDf$site_popularity > 1, ]

# Add TCV column
analysisDf$TCV <- analysisDf$GMV + analysisDf$WMV
phenoCols <- append(regionCols, "TCV")

# Filter for 1 scan per subject
# Sort the dataframe by patient_id and session
analysisDf <- analysisDf[ with(analysisDf, order(analysisDf$scan_id, analysisDf$session)), ]
# Drop all but the first occurrence of each patient_id
analysisDf <- analysisDf[!duplicated(analysisDf$scan_id), ]

names(analysisDf)[names(analysisDf) == "Ventricles"] = "CSF"
row.names(analysisDf) <- NULL


# Save the filtered dataframe
write.csv(analysisDf, paste0(fnBase, "/01_lifespan_subset_filtered.csv"))

# Prep for Combat
# Identify the scan_id + phenotypes to combat
toCombat <- analysisDf[, c('scan_id', phenoCols)]

# Identify covariates to harmonize on/protect
covars <- data.frame(SITE = as.character(analysisDf$site), # <-- harmonize on the first column
                     age_days = analysisDf$age_days, # <-- protect all following columns
                     fs_version = analysisDf$fs_version,
                     sex = analysisDf$sex)

# Save the data and covars dataframes to csvs
write.csv(toCombat, paste0(fnBase, "02_toCombat_phenotypes.csv"), row.names = FALSE)
write.csv(covars, paste0(fnBase, "02_toCombat_covariates.csv"), row.names = FALSE)

# STOP HERE AND RUN COMBAT VIA PYTHON
combatCommand <- paste0("python /Users/youngjm/Projects/mpr_analysis/runNeuroHarmonizeLifespanSubset.py -p ",
                        fnBase, "02_toCombat_phenotypes.csv -c ",
                        fnBase, "02_toCombat_covariates.csv -o ",
                        fnBase, "03_postCombat.csv")
system(combatCommand)

combattedDf <- read.csv(paste0(fnBase, "03_postCombat.csv"))
head(combattedDf)
# combattedDf$scan_id <- as.factor(combattedDf$scan_id)
# analysisDf$scan_id <- as.factor(analysisDf$scan_id)
print(dim(analysisDf))
print(dim(combattedDf))
combattedDf <- merge(combattedDf, analysisDf[, metaCols], by='scan_id')

# Keep only scans with all of the data
combattedDf <- combattedDf[complete.cases(combattedDf), ]
print(dim(combattedDf))

# Save the resulting dataframe with combatted data and metadata
write.csv(combattedDf, paste0(fnBase, '06_combatted_lifespan_subset.csv'), row.names = FALSE)


## ----------------------------------------------------
# Part 3: use GAMLSS to get median curves

# store median centiles somehow here
# Log age column
# combattedDf$site_popularity <- analysisDf$site_popularity
combattedDf <- analysisDf # the data should have already been combatted
combattedDf$logAge <- log(combattedDf$age_days, base=10)
names(combattedDf)[names(combattedDf) == "new_euler"] = "SurfaceHoles"
ageLimited <- sort(read.csv(paste0(fnBase, "age_range.csv"))$x)
# fs_newdata <- case_when(
#   ageLimited <= log(365.25*3+280, base=10) ~ 'FSInfant',
#   TRUE ~ "FS6_T1"
# )
fs_newdata <- c(rep(as.factor("FS6_T1"), length(ageLimited)))
fs_newdata <- as.factor(fs_newdata)



# Set up a list of tick marks to use on log(post-conception age) x-axes
tickMarks <- c()
for (year in c(0, 1, 2, 5, 10, 20)){ # years
  tickMarks <- append(tickMarks, log(year*365.25 + 280, base=10))
}
tickLabels <- c("Birth", "1", "2", "5", "10", "20")

age <- c()
feature <- c()
value <- c()

for (p in phenoCols){
  print(p)
  # For plot 2:
  # 1. Generate GAMLSS models
  formula <- as.formula(paste0(p, "~fp(logAge, npoly=3) + sex - 1")) #+ fs_version
  gamModel <- gamlss(formula = formula,
                    sigma.formula = formula,
                    nu.formula = as.formula(paste0(p, "~1")),
                    family = GG,
                    data = na.omit(combattedDf),
                    control = gamlss.control(n.cyc = 200),  # lifespan
                    trace = F)
  print("finished training the models")
  
  # 2. Predict phenotype values for set age range
  newDataM <- data.frame(logAge=sort(ageLimited),
                         # fs_version=fs_newdata,
                         sex=c(rep(as.factor("Male"),  length(ageLimited))))

  predModelM <- predictAll(gamModel, newdata=newDataM)
  
  newDataF <- data.frame(logAge=sort(ageLimited),
                         # fs_version=fs_newdata,
                         sex=c(rep(as.factor("Female"),  length(ageLimited))))

  predModelF <- predictAll(gamModel, newdata=newDataF)
  
  # Predict the median centile for the M and F models
  phenoMedianPredsM <- qGG(c(0.5), 
                           mu=predModelM$mu, 
                           sigma=predModelM$sigma, 
                           nu=predModelM$nu)
  
  phenoMedianPredsF <- qGG(c(0.5), 
                           mu=predModelF$mu, 
                           sigma=predModelF$sigma, 
                           nu=predModelF$nu)
  phenoMedianPreds <- (phenoMedianPredsF + phenoMedianPredsM)/2
  
  fanCentiles <- c()
  desiredCentiles <- c(0.004, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.996)
  for (i in c(1:length(desiredCentiles))){
    print(desiredCentiles[[i]])
    print(i)
    phenoMedianPredsM <- qGG(desiredCentiles[[i]], 
                             mu=predModelM$mu, 
                             sigma=predModelM$sigma, 
                             nu=predModelM$nu)
    
    phenoMedianPredsF <- qGG(desiredCentiles[[i]], 
                             mu=predModelF$mu, 
                             sigma=predModelF$sigma, 
                             nu=predModelF$nu)
    fanCentiles[[i]] <- (phenoMedianPredsF + phenoMedianPredsM)/2
  }
  
  # 5. Plot CLIP vs Lifespan
  plot <- ggplot() +
    geom_point(data=combattedDf, aes(x=logAge, y=combattedDf[,p], color=sex), alpha=0.3) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[1]]), alpha=0.2) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[2]]), alpha=0.4) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[3]]), alpha=0.6) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[4]]), alpha=0.8) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[5]])) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[6]]), alpha=0.8) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[7]]), alpha=0.6) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[8]]), alpha=0.4) +
    geom_line(aes(x=sort(ageLimited), y=fanCentiles[[9]]), alpha=0.2) +
    # scale_linetype_manual(values = c('solid'), name="50th Centile")+
    scale_x_continuous(breaks=tickMarks, labels=tickLabels,
                       limits=c(tickMarks[[1]], max(ageLimited))) +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(subtitle = p) +
    xlab("Age at scan (log(years))") +
    ylab("Phenotype Value") + 
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 18))
  
  png(file=paste0(fnBase, "fig_growth_chart_", p, ".png"),
      width=600, height=400)
  print(plot)
  dev.off()
  
  # Prepare values to save as csv
  age <- append(age, ageLimited)
  feature <- append(feature, rep(p, length(ageLimited)))
  value <- append(value, phenoMedianPreds)
}

slipAgeMedianDf <- data.frame(age=age, feature=feature, value=value)
write.csv(slipAgeMedianDf, paste0(fnBase, "slip-agemedian_centiles.csv"), row.names = FALSE)

print(table(combattedDf$study))
print(table(combattedDf[combattedDf$sex == "Female", ]$study))
print(table(combattedDf[combattedDf$sex == "Male", ]$study))
print(table(combattedDf$sex))

combattedDf$study <- as.factor(combattedDf$study)

for (s in levels(combattedDf$study)) {
  print(s)
  print((min(combattedDf[combattedDf$study == s,]$age_days - 280))/365.25)
  print((max(combattedDf[combattedDf$study == s,]$age_days - 280))/365.25)
  print((median(combattedDf[combattedDf$study == s,]$age_days - 280))/365.25)
}

print((min(combattedDf$age_days - 280))/365.25)
print((max(combattedDf$age_days - 280))/365.25)
print((median(combattedDf$age_days - 280))/365.25)

