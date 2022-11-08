gc()
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")
library(dplyr)
library(ggseg)
library(gridExtra)

library(ggplot2)
library(stringr)
library(patchwork) # graph organization within a figure

library(gamlss) #to fit model
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")

## Step 7: Prepare regional phenotype data -----------------------

# Load the data
# fn <- '/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
# fn <- '/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
fn <- '/Users/youngjm/Data/clip/fs6_stats/06_combatted_fs_plus_metadata.csv'
t <- "Age at Peak of Original Regional Phenotypes"
fnOut <- '/Users/youngjm/Data/clip/figures/2022-10-18_age_at_peak.png'

brainDf <- read.csv(fn)
brainDf$sex <- as.factor(brainDf$sex)
brainDf$scanner_id <- as.factor(brainDf$scanner_id)
brainDf$scan_reason_categories <- as.factor(brainDf$scan_reason_categories)
brainDfCols <- colnames(brainDf)

regionalPhenotypes <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal',
                        # 'corpus callosum', 
                        'cuneus', # 'entorhinal', 'frontal pole', # these two were excluded from Lifespan
                        'fusiform', 'inferior parietal', 'inferior temporal', 'insula',
                        'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal',
                        'lingual', 'medial orbitofrontal', 'middle temporal', 'paracentral',
                        'parahippocampal', 'pars opercularis', 'pars orbitalis',
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
                              'parahippocampal', 'parsopercularis', 'parsorbitalis',
                              'parstriangularis', 'pericalcarine', 'postcentral',
                              'posteriorcingulate', 'precentral', 'precuneus', 
                              'rostralanteriorcingulate', 'rostralmiddlefrontal',
                              'superiorfrontal', 'superiorparietal', 'superiortemporal',
                              'supramarginal', # 'temporalpole', 
                              'transversetemporal')


grayVolCols <- brainDfCols[endsWith(brainDfCols, '_grayVol')]
print(grayVolCols)
phenoCols <- c()
for (pheno in parsedRegionalPhenotypes) {
  cols <- grayVolCols[grepl(paste("_",pheno, sep=''), grayVolCols, fixed=TRUE)]
  phenoCols <- c(phenoCols, cols)
}
grayVolLhCols <- sort(phenoCols[startsWith(phenoCols, 'lh_')])
grayVolRhCols <- sort(phenoCols[startsWith(phenoCols, 'rh_')])

# Let's calculate the bilateral average of each phenotype
for (i in 1:length(parsedRegionalPhenotypes)){
  brainDf[[parsedRegionalPhenotypes[[i]]]] <- (brainDf[[grayVolLhCols[[i]]]] + brainDf[[grayVolRhCols[[i]]]])/2
}

globalPhenoCols <- c('TotalGrayVol', 'CerebralWhiteMatterVol', 'SubCortGrayVol',
                     'eTIV', 'VentricleVolume', 'CorticalSurfaceArea', 
                     'MeanCorticalThickness', 'TCV')
regionalMetaCols <- setdiff(brainDfCols, c(globalPhenoCols, grayVolCols))
grayVolDf <- brainDf[append(regionalMetaCols, parsedRegionalPhenotypes)]
grayVolDf$logAge <- log(grayVolDf$age_at_scan_days+280, 10)


## Results 4: Make regional phenotype plots for age at peak --------------------
##
# Make the composite plot and tables for a given dataframe/set of phenotypes
# @param df A dataframe where rows are scans and columns are phenotypes and features
# @param phenotypes A vector of phenotypes to analyze
# @param colNames A vector of cleaned names to use for the table columns
# @param rowNames A vector of cleaned names to use for the table rows
# @param title A string to differentiate the figures for this dataframe vs. other dataframes
# generateHemispherePlots <- function(df, phenotypes, phenotypesNoWS, title, sex=""){
  # Initialize variables
  # Set up empty variables to generate table later

ageAtPeak <- c()

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
  
  newDataM <- data.frame(logAge=sort(ageLimited),
                         SurfaceHoles=c(rep(median(grayVolDf$SurfaceHoles), length(ageLimited))),
                         sex=c(rep(as.factor("M"),  length(ageLimited))))
  clipPredModelM <- predictAll(gamModel, newdata=newDataM)
  
  newDataF <- data.frame(logAge=sort(ageLimited),
                         SurfaceHoles=c(rep(median(grayVolDf$SurfaceHoles), length(ageLimited))),
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

  # Get the age at peak - THE PROBLEM IS HERE
  # print(clipPredModel)
  
  maxIdx <- which.max(phenoMedianPreds)
  # Convert ageLimited[maxIdx] from log(age_in_days+280) to years
  peakAge <- (10^(ageLimited[maxIdx])-280)/365.25
  ageAtPeak <- append(ageAtPeak, c(peakAge))
  
  plt <- ggplot() + 
    geom_point(aes(x=logAge, y=grayVolDf[, phenotype]), data=grayVolDf, color="red") + 
    geom_smooth(aes(x=ageLimited, y=phenoMedianPreds), color="blue") +
    geom_vline(xintercept=ageLimited[maxIdx], color="blue") +
    geom_hline(yintercept=max(phenoMedianPreds)) +
    # geom_point(aes(x=ageLimited, y=phenoMedianPreds), color="blue") +
    labs(title=paste(phenotype, "Max volume:", max(phenoMedianPreds), "at age", peakAge)) + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
    # geom_vline(xintercept=)
  print(plt)
}

# Convert the lists into a dataframe
results = as.data.frame(cbind(region=regionalPhenotypes, 
                              feat=parsedRegionalPhenotypes,
                              em=as.numeric(ageAtPeak)),
                        stringsAsFactors=F)

print(results)

# Plot the dataframe using ggseg
p <- results %>% 
  ggseg(mapping=aes(fill=as.numeric(em)),
        hemisphere='left') +
  labs(title = "Age at Peak (years)") +
  theme(axis.title = element_blank()) +
  scale_fill_gradient(low = "blue", high = "red", na.value = NA, name="Age") 

grid.arrange(p)

regionalPeaks <- results #generateHemispherePlots(grayVolDf, regionalPhenotypes, parsedRegionalPhenotypes, '')

# LOH
lifespanDf <- read.csv('/Users/youngjm/Data/lifespan_growth_charts/Lifespan_Data_Peaks_Table_2_2.csv')
comboDf <- merge(regionalPeaks, lifespanDf, by='feat')
comboDf <- comboDf[comboDf$feat %in% parsedRegionalPhenotypes, ]
comboDf$clipPeak <- as.numeric(comboDf$em)
comboDf$peakDiff <- as.numeric(comboDf$Peak) - as.numeric(comboDf$clipPeak)
minAge <- min(comboDf$clipPeak, comboDf$Peak)
maxAge <- max(comboDf$clipPeak, comboDf$Peak)
r <- cor(comboDf$clipPeak, comboDf$Peak)

avgVol <- c()
for (pheno in comboDf$feat){
  print(pheno)
  avg <- mean(grayVolDf[[pheno]])/4000
  avgVol <- append(avgVol, avg)
}
comboDf$avgVol <- avgVol

plots <- c()
plots[[1]] <- comboDf %>%
  ggseg(mapping=aes(fill=clipPeak),
        hemisphere='left') +
  labs(title = "CLIP") +
  scale_fill_gradient(low = "blue", high = "red", 
                      limits=c(minAge, maxAge),
                      na.value = NA, name="Age (years)") 

plots[[2]] <- comboDf %>%
  ggseg(mapping=aes(fill=Peak),
        hemisphere='left') +
  labs(title = "Lifespan") +
  scale_fill_gradient(low = "blue", high = "red", 
                      limits=c(minAge, maxAge),
                      na.value = NA, name="Age (years)") 

comboDf$logClipPeak <- log(comboDf$clipPeak*365.25+280, base=10)
comboDf$logLifespanPeak <- log(comboDf$Peak*365.25+280, base=10)
plots[[3]] <- ggplot(data=comboDf, aes(color=as.factor(region), shape=as.factor(region), fill=as.factor(region))) +
  geom_point(mapping = aes(x=logClipPeak, y=logLifespanPeak), size=avgVol) +
  geom_abline(slope = 1) +
  scale_shape_manual(values = rep(21:25, 7), name="Region") +
  scale_color_manual(values = rep(cbbPalette, 4), name="Region") +
  scale_fill_manual(values = rep(cbbPalette, 4), name="Region") +
  theme_minimal() +
  # guides(fill = 'none') +
  ylab('Lifespan Age at Peak') + 
  xlab('CLIP Age at Peak') +
  xlim(min(comboDf$logClipPeak, comboDf$logLifespanPeak), max(comboDf$logClipPeak, comboDf$logLifespanPeak)) +
  ylim(min(comboDf$logClipPeak, comboDf$logLifespanPeak), max(comboDf$logClipPeak, comboDf$logLifespanPeak)) +
  labs(title = paste0('Lifespan vs. CLIP (r=', format(r, digits=4), ')')) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

layout <-"
A
B
B
"
brainPlots <- wrap_plots(plots[[1]] + plots[[2]], guides = "collect")
patch <- wrap_plots(brainPlots + plots[[3]] + plot_layout(design=layout), guides="auto")
png(file=fnOut,
    width=750, height=600)
print(patch + plot_annotation(title="Lifespan and CLIP Age at Peak Region Volume"))
dev.off()
