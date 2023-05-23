gc()
# library(dplyr)
library(ggseg)
library(gridExtra)

# library(ggplot2)
library(stringr)
library(patchwork) # graph organization within a figure
library(ggthemes)

# library(gamlss) #to fit model

setwd("/Users/youngjm/Projects/mpr_analysis/r/")
source("lib_mpr_analysis.r")


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",  
                "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

#-------------------------------------------------------------------------------
# Filepaths to change
#-------------------------------------------------------------------------------
fn <- '/Users/youngjm/Data/slip/fs6_stats/06_combatted_fs_plus_metadata.csv'
fnLifespan <- '/Users/youngjm/Data/lifespan_growth_charts/Lifespan_Data_Peaks_Table_2_2.csv'
t <- "Age at Peak of Original Regional Phenotypes"
fnOut <- '/Users/youngjm/Data/slip/figures/2022-10-18_age_at_peak.png'
fnOutSlipRegionalAgeAtPeak <- "/Users/youngjm/Data/slip/fs6_stats/2022-11-21_age_at_peak_slip_lifespan_data.csv"

#-------------------------------------------------------------------------------
# Load the data
#-------------------------------------------------------------------------------
brainDf <- read.csv(fn)
brainDf$sex <- as.factor(brainDf$sex)
brainDf$scanner_id <- as.factor(brainDf$scanner_id)
brainDf$scan_reason_categories <- as.factor(brainDf$scan_reason_categories)
brainDfCols <- colnames(brainDf)

#-------------------------------------------------------------------------------
# Extract the regional phenotype data for the left and right hemispheres
#-------------------------------------------------------------------------------
# Set up lists of regional phenotypes
regionalPhenotypes <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal',
                        # 'corpus callosum', 
                        'cuneus', # 'entorhinal', 'frontal pole', # these two were excluded from Lifespan
                        'fusiform', 'inferior parietal', 'inferior temporal', 'insula',
                        'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal',
                        'lingual', 'medial orbitofrontal', 'middle temporal', 'paracentral',
                        #'parahippocampal', # did not peak in the SLIP age range
                        'pars opercularis', 'pars orbitalis',
                        'pars triangularis', 'pericalcarine', 'postcentral',
                        'posterior cingulate', 'precentral', 'precuneus', 
                        'rostral anterior cingulate', 'rostral middle frontal',
                        'superior frontal', 'superior parietal', 'superior temporal',
                        'supramarginal', #'temporal pole', # excluded from Lifespan
                        'transverse temporal')
parsedRegionalPhenotypes <- c('bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal',
                              # 'corpus callosum', 
                              'cuneus', #'entorhinal', 'frontalpole', 
                              'fusiform', 'inferiorparietal', 'inferiortemporal', 'insula',
                              'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal',
                              'lingual', 'medialorbitofrontal', 'middletemporal', 'paracentral',
                              #'parahippocampal',
                              'parsopercularis', 'parsorbitalis',
                              'parstriangularis', 'pericalcarine', 'postcentral',
                              'posteriorcingulate', 'precentral', 'precuneus', 
                              'rostralanteriorcingulate', 'rostralmiddlefrontal',
                              'superiorfrontal', 'superiorparietal', 'superiortemporal',
                              'supramarginal', # 'temporalpole', 
                              'transversetemporal')

# Identify the regional volumetric phenotype columns
grayVolCols <- brainDfCols[endsWith(brainDfCols, '_grayVol')]
# print(grayVolCols)
phenoCols <- c()
for (pheno in parsedRegionalPhenotypes) {
  cols <- grayVolCols[grepl(paste("_",pheno, sep=''), grayVolCols, fixed=TRUE)]
  phenoCols <- c(phenoCols, cols)
}
grayVolLhCols <- sort(phenoCols[startsWith(phenoCols, 'lh_')])
grayVolRhCols <- sort(phenoCols[startsWith(phenoCols, 'rh_')])

# Calculate the bilateral average of each phenotype
for (i in 1:length(parsedRegionalPhenotypes)){
  brainDf[[parsedRegionalPhenotypes[[i]]]] <- (brainDf[[grayVolLhCols[[i]]]] + brainDf[[grayVolRhCols[[i]]]])/2
}


# Pull only the metadata columns and the bilateral averages of each phenotype
globalPhenoCols <- c('TotalGrayVol', 'CerebralWhiteMatterVol', 'SubCortGrayVol',
                     'eTIV', 'VentricleVolume', 'CorticalSurfaceArea', 
                     'MeanCorticalThickness', 'TCV')
regionalMetaCols <- setdiff(brainDfCols, c(globalPhenoCols, grayVolCols))
grayVolDf <- brainDf[append(regionalMetaCols, parsedRegionalPhenotypes)]
grayVolDf$logAge <- log(grayVolDf$age_at_scan_days+280, 10)


#-------------------------------------------------------------------------------
# Calculate the age at peak volume for each region
#-------------------------------------------------------------------------------
ageAtPeak <- c()

# Set the x axis limits
minLogAge <- min(grayVolDf$logAge)
maxLogAge <- max(grayVolDf$logAge)
ageLimited <- seq(from=minLogAge, to=maxLogAge, by=0.005)

# for phenotype in phenotypes...
for (phenotype in parsedRegionalPhenotypes){
  print(phenotype)
  form <- as.formula(paste0(phenotype, "~fp(logAge, npoly=3) + SurfaceHoles + sex -1")) 
  gamModel <- gamlss(formula = form, 
                    sigma.formula = form,
                    nu.formula = as.formula(paste0(phenotype, "~1")),
                    family = GG, 
                    data=na.omit(grayVolDf), 
                    control = gamlss.control(n.cyc = 200),  # lifespan
                    trace = F)
  
  phenoMedianPreds <- predictCentilesForAgeRange(gamModel, ageLimited, 
                                                 euler=median(grayVolDf$SurfaceHoles))

  # Get the age at peak 
  maxIdx <- which.max(phenoMedianPreds)
  # Convert ageLimited[maxIdx] from log(age_in_days+280) to years
  peakAge <- (10^(ageLimited[maxIdx])-280)/365.25
  ageAtPeak <- append(ageAtPeak, c(peakAge))
  
  # User can remove this section if they would like - it shows age at peak for 
  # for each phenotype on the growth chart for that phenotype
  plt <- ggplot() + 
    geom_point(aes(x=logAge, y=grayVolDf[, phenotype]), data=grayVolDf, color="red") + 
    geom_smooth(aes(x=ageLimited, y=phenoMedianPreds), color="blue") +
    geom_vline(xintercept=ageLimited[maxIdx], color="blue") +
    geom_hline(yintercept=max(phenoMedianPreds)) +
    labs(title=paste(phenotype, "Max volume:", max(phenoMedianPreds), "at age", peakAge)) + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  print(plt)
}

# Plot the age at peak volume in the regional brain graph (ggseg)
# Convert the lists into a dataframe
results = as.data.frame(cbind(region=regionalPhenotypes, 
                              feat=parsedRegionalPhenotypes,
                              em=as.numeric(ageAtPeak)),
                        stringsAsFactors=F)

print(results)

# Plot the dataframe using ggseg
p <- results %>% 
  ggseg(mapping=aes(fill=as.numeric(em))) +
  labs(title = "Age at Peak (years)") +
  theme(axis.title = element_blank()) +
  scale_fill_gradient_tableau(palette = "Blue", na.value = NA, name="Age") 

grid.arrange(p)

regionalPeaks <- results

#-------------------------------------------------------------------------------
# Build the composite figure showing slip and lifespan age at peak volume for 
# regional phenotypes 
#-------------------------------------------------------------------------------
lifespanDf <- read.csv(fnLifespan)
comboDf <- merge(regionalPeaks, lifespanDf, by='feat')
comboDf <- comboDf[comboDf$feat %in% parsedRegionalPhenotypes, ]
comboDf$slipPeak <- as.numeric(comboDf$em)
comboDf$peakDiff <- as.numeric(comboDf$Peak) - as.numeric(comboDf$slipPeak)
minAge <- min(comboDf$slipPeak, comboDf$Peak)
maxAge <- max(comboDf$slipPeak, comboDf$Peak)
r <- cor(comboDf$slipPeak, comboDf$Peak)

avgVol <- c()
for (pheno in comboDf$feat){
  print(pheno)
  avg <- mean(grayVolDf[[pheno]])/4000
  avgVol <- append(avgVol, avg)
}
comboDf$avgVol <- avgVol

write.csv(comboDf, fnOutSlipRegionalAgeAtPeak)

# Make the first plot: SLIP age at peak volume
plots <- c()
plots[[1]] <- comboDf %>%
  ggseg(mapping=aes(fill=slipPeak),
        hemisphere='left') +
  labs(title = "SLIP") +
  scale_fill_gradient_tableau(palette = "Blue",
                              limits=c(minAge, maxAge),
                              na.value = NA, name="Age (years)") + #, trans="log") +
  theme_void()
  # theme_brain(text.size=14, text.family="sans")

# Make the second plot: Lifespan age at peak volume
plots[[2]] <- comboDf %>%
  ggseg(mapping=aes(fill=Peak),
        hemisphere='left') +
  labs(title = "LBCC") +
  scale_fill_gradient_tableau(palette = "Blue",
                              limits=c(minAge, maxAge),
                              na.value = NA, name="Age (years)") + #, trans="log") +
  theme_void()
  # theme_brain(text.size=14, text.family="sans")

# Make the third plot: age vs age (correlation)
# Adding tick mark labels for log(post conception age in days)
tickMarks <- c()
for (year in c(2, 5, 10, 20)){ # years
  tickMarks <- append(tickMarks, log(year*365.25 + 280, base=10))
}
tickLabels <- c("2", "5", "10", "20")

comboDf$logSlipPeak <- log(comboDf$slipPeak*365.25+280, base=10)
comboDf$logLifespanPeak <- log(comboDf$Peak*365.25+280, base=10)
plots[[3]] <- ggplot(data=comboDf, aes(color=as.factor(region), shape=as.factor(region), fill=as.factor(region))) +
  geom_point(mapping = aes(x=logSlipPeak, y=logLifespanPeak), size=avgVol) +
  geom_abline(slope = 1) +
  scale_shape_manual(values = rep(21:25, 7), name="Region") +
  scale_color_manual(values = rep(cbbPalette, 5), name="Region") +
  scale_fill_manual(values = rep(cbbPalette, 5), name="Region") +
  theme_minimal() +
  ylab('LBCC Age at Peak (log(years))') + 
  xlab('SLIP Age at Peak (log(years))') +
  scale_x_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(comboDf$logSlipPeak, comboDf$logLifespanPeak))) +
  scale_y_continuous(breaks=tickMarks, labels=tickLabels, 
                     limits=c(tickMarks[[1]], max(comboDf$logSlipPeak, comboDf$logLifespanPeak))) +
  labs(title = paste0("LBCC vs. SLIP (Spearman's r = ", format(r, digits=4), ')')) + 
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size=14),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

layout <-"
A
B
B
"
brainPlots <- wrap_plots(plots[[1]] + plots[[2]], guides = "collect")
patch <- wrap_plots(brainPlots + plots[[3]] + plot_layout(design=layout))
png(file=fnOut,
    width=1000, height=600)
print(patch + plot_annotation(tag_levels = 'A')) 
#+ plot_annotation(title="LBCC and SLIP Age at Peak Region Volume") &
        #theme(plot.title = element_text(size=16)))
dev.off()

