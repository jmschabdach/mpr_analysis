gc()
dev.off(dev.list()["RStudioGD"])

library(ggplot2)
# library(ggpubr)
library(dplyr)
library(mgcv)
library(tidymv)
library(patchwork) # graph organization within a figure
# library(gtsummary)
# library(grid)
# library(harrypotter)
# library(stringr)
# library(gridExtra)
# library(reshape2)
# library(tables)
# library(gridExtra)
# library(data.table)
# library(formattable)
# library(tidyr)
# library(ggseg)

source("lib_mpr_analysis.r")

# Colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")

fn <- '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_combatted_covariates_removed_plus_metadata.csv'
analysisDf <- read.csv(fn)

# Convert some columns to factors
toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 'scanner_id', 'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)


## Step 2: Generate basic demographic plots ------------------------------------
# These plots should be consistent across all analyses

qcColors <- c("#999999", "#0072B2", "#56B4E9")

# Make a plot for the distribution of ages with coloring for sex 
pAgesMale <- generateAgeDistributionPlot(analysisDf, 'M', qcColors)
pAgesFemale <- generateAgeDistributionPlot(analysisDf, 'F', qcColors)

pScannersMale <- generateScannerDistributionPlot(analysisDf, 'M', qcColors)
pScannersFemale <- generateScannerDistributionPlot(analysisDf, 'F', qcColors)

# Make plots for QC
pEulerQcMale <- generateEulerQcDistributionPlot(analysisDf, 'M', qcColors)
pEulerQcFemale <- generateEulerQcDistributionPlot(analysisDf, 'F', qcColors)

pAgeQcMale <- generateAgeQcDistributionPlot(analysisDf, 'M', qcColors)
pAgeQcFemale <- generateAgeQcDistributionPlot(analysisDf, 'F', qcColors)

ttest2v1M <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 1, ]$age_in_years)
ttest1v0M <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 1, ]$age_in_years, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v0M <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v1F <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 1, ]$age_in_years)
ttest1v0F <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 1, ]$age_in_years, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v0F <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v1All <- t.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$rawdata_image_grade == 1, ]$age_in_years)
ttest1v0All <- t.test(analysisDf[analysisDf$rawdata_image_grade == 1, ]$age_in_years, analysisDf[analysisDf$rawdata_image_grade == 0, ]$age_in_years)
ttest2v0All <- t.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$age_in_years, analysisDf[analysisDf$rawdata_image_grade == 0, ]$age_in_years)


qcAgeTableM <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                          pvalues = c(ttest1v0M$p.value, ttest2v0M$p.value, ttest2v1M$p.value),
                          ciLower = c(ttest1v0M$conf.int[1], ttest2v0M$conf.int[1], ttest2v1M$conf.int[1]),
                          ciUpper = c(ttest1v0M$conf.int[2], ttest2v0M$conf.int[2], ttest2v1M$conf.int[2]))

qcAgeTableF <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                          pvalues = c(ttest1v0F$p.value, ttest2v0F$p.value, ttest2v1F$p.value),
                          ciLower = c(ttest1v0F$conf.int[1], ttest2v0F$conf.int[1], ttest2v1F$conf.int[1]),
                          ciUpper = c(ttest1v0F$conf.int[2], ttest2v0F$conf.int[2], ttest2v1F$conf.int[2]))

qcAgeTableAll <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                            pvalues = c(ttest1v0All$p.value, ttest2v0All$p.value, ttest2v1All$p.value),
                            ciLower = c(ttest1v0All$conf.int[1], ttest2v0All$conf.int[1], ttest2v1All$conf.int[1]),
                            ciUpper = c(ttest1v0All$conf.int[2], ttest2v0All$conf.int[2], ttest2v1All$conf.int[2]))

# Test for Euler and QC
ttest2v1EulerM <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles)
ttest1v0EulerM <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v0EulerM <- t.test(analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'M' & analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v1EulerF <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles)
ttest1v0EulerF <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v0EulerF <- t.test(analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$sex == 'F' & analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v1EulerAll <- t.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles)
ttest1v0EulerAll <- t.test(analysisDf[analysisDf$rawdata_image_grade == 1, ]$SurfaceHoles, analysisDf[analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)
ttest2v0EulerAll <- t.test(analysisDf[analysisDf$rawdata_image_grade == 2, ]$SurfaceHoles, analysisDf[analysisDf$rawdata_image_grade == 0, ]$SurfaceHoles)

qcAgeTableEulerM <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                               pvalues = c(ttest1v0EulerM$p.value, ttest2v0EulerM$p.value, ttest2v1EulerM$p.value),
                               ciLower = c(ttest1v0EulerM$conf.int[1], ttest2v0EulerM$conf.int[1], ttest2v1EulerM$conf.int[1]),
                               ciUpper = c(ttest1v0EulerM$conf.int[2], ttest2v0EulerM$conf.int[2], ttest2v1EulerM$conf.int[2]))

qcAgeTableEulerF <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                               pvalues = c(ttest1v0EulerF$p.value, ttest2v0EulerF$p.value, ttest2v1EulerF$p.value),
                               ciLower = c(ttest1v0EulerF$conf.int[1], ttest2v0EulerF$conf.int[1], ttest2v1EulerF$conf.int[1]),
                               ciUpper = c(ttest1v0EulerF$conf.int[2], ttest2v0EulerF$conf.int[2], ttest2v1EulerF$conf.int[2]))

qcAgeTableEulerAll <- data.frame(comparison = c("0 vs. 1", "0 vs. 2", "1 vs. 2"),
                                 pvalues = c(ttest1v0EulerAll$p.value, ttest2v0EulerAll$p.value, ttest2v1EulerAll$p.value),
                                 ciLower = c(ttest1v0EulerAll$conf.int[1], ttest2v0EulerAll$conf.int[1], ttest2v1EulerAll$conf.int[1]),
                                 ciUpper = c(ttest1v0EulerAll$conf.int[2], ttest2v0EulerAll$conf.int[2], ttest2v1EulerAll$conf.int[2]))



# Make a table containing the number of scans with each rating divided by sex
qcCountTable <- data.frame(rating = c("2", "1", "0"),
                           Male = c(dim(analysisDf[analysisDf$sex == "M" & analysisDf$rawdata_image_grade == 2, ])[1],
                                    dim(analysisDf[analysisDf$sex == "M" & analysisDf$rawdata_image_grade == 1, ])[1],
                                    dim(analysisDf[analysisDf$sex == "M" & analysisDf$rawdata_image_grade == 0, ])[1]),
                           Female = c(dim(analysisDf[analysisDf$sex == "F" & analysisDf$rawdata_image_grade == 2, ])[1],
                                      dim(analysisDf[analysisDf$sex == "F" & analysisDf$rawdata_image_grade == 1, ])[1],
                                      dim(analysisDf[analysisDf$sex == "F" & analysisDf$rawdata_image_grade == 0, ])[1]))



grob <- patchworkGrob(pAgesMale + pAgesFemale + 
                        pScannersMale + pScannersFemale +
                        pAgeQcMale + pAgeQcFemale +
                        pEulerQcMale + pEulerQcFemale +
                        plot_layout(guides="collect", ncol = 2))
gridExtra::grid.arrange(grob)

# ^^^ SAVE THIS FIGURE

# Plot the distribution of the primary reasons for a scan
analysisDf %>%
  group_by(scan_reason_primary) %>%
  summarise(scan_count = n()) %>%
  ggplot(aes(x=reorder(scan_reason_primary, (scan_count)), y=scan_count)) +
  geom_bar(stat='identity', show.legend = FALSE) +
  theme_minimal() +
  coord_flip() +
  xlab("Primary Reason for Scan") +
  ylab("Number of Scans") +
  ggtitle("Number of Scans Performed for Reason")


# Plot the distribution of top 5 reasons for a scan
analysisDf <- addPrimaryScanReasonCol(analysisDf)
reasonsTable <- summary(analysisDf$top_scan_reason_factors)
pieLabels <- paste(names(reasonsTable), reasonsTable, sep='\n')
pie(reasonsTable, pieLabels, col=cbbPalette, mai=c(0,0,0,0),
    main="Top Scan Reasons for Clinical Controls\n(with sample sizes)")

