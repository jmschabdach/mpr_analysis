gc()
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)
library(stringr)
library(patchwork) # graph organization within a figure

library(gamlss) #to fit model
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")

# Part 0: load the CLIP data from multiple files and make sure all of the data
#  types match
fn <- '/Users/youngjm/Data/clip/fs6_stats/07_fully_filtered_postcombat_clip_fs.csv'
clipDf <- read.csv(fn)
fn <- '/Users/youngjm/Data/clip/fs6_stats/fs_gamlss_centiles.csv'
violinDf <- read.csv(fn) 

# Prep CLIP FS data
# Add a column to clip for log(post conception age)
clipDf$logAge <- log(clipDf$age_at_scan_days+280, 10)
# Make sex and scanner_id factors
clipDf$sex <- as.factor(clipDf$sex)
# clipDf$scanner_id <- as.factor(clipDf$scanner_id)
clipDf$proc_ord_year <- as.factor(clipDf$proc_ord_year)
clipDf$top_scan_reason_factors <- as.factor(clipDf$top_scan_reason_factors)
# Drop rows with missing values
clipDf <- drop_na(clipDf)

# Test: filtering out any scans processed using IFS
clipDf <- clipDf[clipDf$age_at_scan_days >= (365.25*3), ]

# Load the Lifespan models
lifespanFnBase <- '/Users/youngjm/Data/clip/fs6_stats/lifespan_'
phenoCaps <- c("GMV", "WMV", "sGMV", "CSF", "TCV")
commonCols <- c("participant", "Age", "age_days", "sex", "study", "fs_version", 
                # "scanner_id", 
                "top_scan_reason_factors", "SurfaceHoles", 
                "GMV", "WMV", "sGMV", "Ventricles", "SA", "CT", "TCV")
phenoColNames <- c()
phenoColValues <- c()
for (p in phenoCaps){
  fn <- paste0(lifespanFnBase, p, ".csv")
  tmpDf <- read.csv(fn)
  if (p == "CSF") { p <- "Ventricles" }
  adjustedPhenoCol <- paste0(p, "Transformed.normalised")
  phenoCentileCol <- paste0(p, "Transformed.q.wre")
  phenoColNames <- append(phenoColNames, c(adjustedPhenoCol, phenoCentileCol))
  phenoColValues[[adjustedPhenoCol]] <- tmpDf[, adjustedPhenoCol]
  phenoColValues[[phenoCentileCol]] <- tmpDf[, phenoCentileCol]
}
lifespanDf <- tmpDf[, commonCols]
lifespanDf[phenoColNames] <- phenoColValues
names(lifespanDf)[names(lifespanDf) == "VentriclesTransformed.normalized"] <- "CSFTransformed.normalized"
names(lifespanDf)[names(lifespanDf) == "VentriclesTransformed.q.wre"] <- "CSFTransformed.q.wre"
names(lifespanDf)[names(lifespanDf) == "age_days"] <- "age_at_scan_days"


# Prep Lifespan data
# Formatting for sex
lifespanDf$sex <- as.factor(lifespanDf$sex)
# lifespanDf$scanner_id <- as.factor(lifespanDf$scanner_id)
# lifespanDf$proc_ord_year <- as.factor(lifespanDf$proc_ord_year)
lifespanDf$top_scan_reason_factors <- as.factor(lifespanDf$top_scan_reason_factors)
# # Scale the Lifespan values - they're N x 10,000
# lc$value <- lc$value * 10000

# Filter the lifespanDf to only include data in the clipDf
lifespanDf <- subset(lifespanDf, participant %in% clipDf$scan_id)
violinDf <- subset(violinDf, vScanId %in% clipDf$scan_id)

# Set up the lifespanViolinDf
regions <- c()
centilesList <- c()
reasons <- c()
participant <- c()
for (p in phenoCaps){
  tmp <- lifespanDf[, grepl(paste0("^", p, "Transformed.q.wre"), names(lifespanDf))]
  centilesList <- append(centilesList, tmp)
  regions <- append(regions, rep(p, dim(lifespanDf)[1]))
  reasons <- append(reasons, lifespanDf$top_scan_reason_factors)
  participant <- append(participant, lifespanDf$participant)
}
centilesLifespanDf <- data.frame(participant, regions, centilesList, reasons)

# Set up variables needed for the analysis
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#D55E00", "#CC79A7", "#F0E442", "#0072B2")
reasonLabels <- c("Developmental disorder", "Clinical eye or vision finding", "Headache", "Other", "Suspected seizure")


# Make the violin plots
violin <- c()
violin[[1]] <- ggplot(data=violinDf, aes(regions, centilesList)) +
  geom_violin(color="gray", fill="gray", alpha=0.5) +
  geom_jitter(height = 0, width=0.15, aes(color=reasons), alpha=0.65) +
  scale_color_manual(values = cbbPalette, labels=reasonLabels, name = "Reason for Scan") +
  labs(title="Distributions of FS Centiles (GAMLSS)") + 
  xlab("Tissue Type") +
  ylab("Centile") + 
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))


violin[[2]] <- ggplot(data=centilesLifespanDf, aes(regions, centilesList)) +
  geom_violin(color="gray", fill="gray", alpha=0.5) +
  geom_jitter(height = 0, width=0.15, aes(color=reasons), alpha=0.65) +
  scale_color_manual(values = cbbPalette, labels=reasonLabels, name = "Reason for Scan") +
  labs(title="Distributions of FS Centiles (Lifespan Models)") + 
  xlab("Tissue Type") +
  ylab("Centile") + 
  theme(axis.line = element_line(colour = "black"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18))

vsPlots <- c()
for (p in c(1:length(phenoCaps))){
  print(p)
  violinDfSmol <- violinDf[regions == phenoCaps[[p]], ]
  lifespanDfSmol <- centilesLifespanDf[regions == phenoCaps[[p]], ]
  violinDfSmol <- violinDfSmol[order(violinDfSmol$vScanId), ]
  lifespanDfSmol <- lifespanDfSmol[order(lifespanDfSmol$participant), ]
  r <- cor(violinDfSmol$centilesList, lifespanDfSmol$centilesList)

  vsPlots[[p]] <- ggplot() +
    geom_point(aes(x=violinDfSmol$centilesList, y=lifespanDfSmol$centilesList, color=lifespanDf$top_scan_reason_factors), alpha=0.3) +
    scale_color_manual(values = cbbPalette, labels=reasonLabels, name = "Reason for Scan") +
    labs(subtitle = paste0("GAMLSS vs. LBCC for ", phenoCaps[[p]], " (r=", format(r, digits=3), ")")) +
    xlab("GAMLSS Centiles") +
    ylab("LBCC Centiles") + 
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 18))
  png(file=paste0("/Users/youngjm/Data/clip/figures/2022-12-12_gamlss_lifespan_", phenoCaps[[p]], ".png"),
                  width=600, height=400)
  print(vsPlots[[p]])
  dev.off()

}
# lower <- wrap_plots(vsPlots, guides="collect")
# print(lower)

top <- wrap_plots(violin, nrow=2, guides = "collect")
png(file="/Users/youngjm/Data/clip/figures/2022-12-12_gamlss_lifespan_centiles.png",
    width=1200, height=800)
print(top)
dev.off()

