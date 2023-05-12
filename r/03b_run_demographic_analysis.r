gc()

library(ggplot2)
library(dplyr)
library(mgcv)
library(stringr)
library(tidymv)
library(patchwork) # graph organization within a figure

source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")

# Colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")

fn <- '/Users/youngjm/Data/slip/fs6_stats/original_phenotypes_demographics.csv'
# fn <- '/Users/youngjm/Data/slip/fs6_stats/06_combatted_plus_metadata.csv'
analysisDf <- read.csv(fn)
fnOut <- '/Users/youngjm/Data/slip/figures/2022-11-07_demographics_figure.png'

# Convert some columns to factors
toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 'scanner_id', 'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

# names(analysisDf)[names(analysisDf) == "TotalGrayVol"] <- 'GMV'
# names(analysisDf)[names(analysisDf) == "CerebralWhiteMatterVol"] <- 'WMV'
# names(analysisDf)[names(analysisDf) == "SubCortGrayVol"] <- 'sGMV'
# names(analysisDf)[names(analysisDf) == "VentricleVolume"] <- 'CSF'
# names(analysisDf)[names(analysisDf) == "CorticalSurfaceArea"] <- 'SA'
# names(analysisDf)[names(analysisDf) == "MeanCorticalThickness"] <- 'CT'

demoFigs <- c()

analysisDf$scanner_id <-gsub("Scanner ","",as.character(analysisDf$scanner_id))
# freqTable <- as.data.frame(table(analysisDf$scanner_id))
# analysisDf$scannerId <- factor(analysisDf$scanner_id, 
#                                levels = freqTable$Var1[order(freqTable$Freq, decreasing = T)])


# Make a bar chart for age at scan
demoFigs[[1]] <- ggplot(analysisDf, aes(x=age_in_years, fill=sex)) +
  # geom_bar(aes(x=age_in_years, fill=sex), width=0.25, position="stack") +
  geom_histogram(position=position_stack(), binwidth = 0.5) +
  scale_fill_manual(values = cbbPalette[2:3], name = "Sex") +
  labs(title="Distribution of Age at Scan") + 
  ylab("Number of Scans") +
  xlab("Age at Scan (years)") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14)) 

# Make a bar chart for age at scan
scannerLabels <- str_pad(c(1:length(levels(analysisDf$scanner_id))), 2, pad = "0")
scannerLabels <- paste('Scanner', scannerLabels, sep='\n')

demoFigs[[2]] <- ggplot(analysisDf, aes(x=scanner_id, fill=sex)) +
  geom_bar(position="stack") +
  # geom_histogram(position=position_stack(), binwidth = 0.5) +
  scale_fill_manual(values = cbbPalette[2:3], name = "Sex") +
  labs(title="Scanner Distributions") + 
  ylab("Number of Scans") +
  xlab("Scanner ID") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14)) 


# Plot age vs. average image quality for male and female 
test <- aov(average_grade~age_in_years, data=analysisDf)
pval <- summary(test)[[1]][["Pr(>F)"]][1]
fval <- summary(test)[[1]][["F value"]][1]
dof <- summary(test)[[1]][["Df"]][1]
r <- cor(analysisDf$average_grade, analysisDf$age_in_years)
demoFigs[[3]] <- ggplot(analysisDf) +
  geom_hline(yintercept=1.0, color="black", alpha=0.5) +
  # geom_point(alpha=0.75, aes(x=age_in_years, y=average_grade, color=sex)) +
  geom_jitter(height = 0.025, width=0.01, alpha=0.5, aes(x=age_in_years, y=average_grade, color=sex)) +
  scale_color_manual(values = cbbPalette[2:3], name = "Sex") +
  labs(title="Average QC Rating vs. Age at Scan",
       subtitle = paste0("Pearson's r = ", format(r, digits=3))) +
       # subtitle = paste0(" (F(",paste(dim(analysisDf)[1]-1), ", ", paste(dof),") = ",
       #              format(fval, digits=3), ", p < ", format(pval, digits=3), ")")) + 
  ylab("Average QC Grade") +
  xlab("Age at Scan (years)") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="gray"), # element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14)) 

# Plot scanner ID vs. average image quality for male and female 
demoFigs[[4]] <- ggplot(analysisDf) +
  geom_hline(yintercept=1.0, color="black", alpha=0.5) +
  # geom_point(alpha=0.75, aes(x=age_in_years, y=average_grade, color=sex)) +
  geom_jitter(height = 0.025, width=0.25, alpha=0.5, aes(x=scanner_id, y=average_grade, color=sex)) +
  scale_color_manual(values = cbbPalette[2:3], name = "Sex") +
  labs(title="Average QC Rating vs. Scanner ID") + 
  ylab("Average QC Grade") +
  xlab("Scanner ID") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="gray"), # element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14)) 

# Plot scanner ID vs. Euler for male and female 
test <- aov(average_grade~SurfaceHoles, data=analysisDf)
pval <- summary(test)[[1]][["Pr(>F)"]][1]
fval <- summary(test)[[1]][["F value"]][1]
dof <- summary(test)[[1]][["Df"]][1]
r <- cor(analysisDf$average_grade, analysisDf$SurfaceHoles)
demoFigs[[5]] <- ggplot(analysisDf) +
  geom_hline(yintercept=1.0, color="black", alpha=0.5) +
  # geom_point(alpha=0.75, aes(x=SurfaceHoles, y=average_grade, color=sex)) +
  geom_jitter(height = 0.1, width=0.1, alpha=0.5, aes(x=SurfaceHoles, y=average_grade, color=sex)) +
  scale_color_manual(values = cbbPalette[2:3], name = "Sex") +
  labs(title=paste0("Average QC Rating vs. Euler Number"),
       subtitle = paste0("Pearson's r = ", format(r, digits=3))) + 
       # subtitle = paste0(" (F(",paste(dim(analysisDf)[1]-1), ", ", paste(dof),") = ",
       #                  format(fval, digits=3), ", p < ", format(pval, digits=3), ")")) + 
  ylab("Average QC Grade") +
  xlab("Euler Number") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="gray"), # element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14)) 

# Plot Euler vs. average image quality for male and female
demoFigs[[6]] <- ggplot(analysisDf) +
  # geom_point(alpha=0.75, aes(x=age_in_years, y=average_grade, color=sex)) +
  geom_jitter(height = 0.025, width=0.25, alpha=0.5, aes(x=scanner_id, y=SurfaceHoles, color=sex)) +
  scale_color_manual(values = cbbPalette[2:3], name = "Sex") +
  labs(title="Euler Number vs. Scanner ID") + 
  ylab("Euler Number") +
  xlab("Scanner ID") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="gray"), # element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 14)) 

# # Plot age vs. Euler for male and female 
# ageEulerSex <- ggplot(analysisDf) +
#   geom_hline(yintercept=1.0, color="black", alpha=0.5) +
#   geom_point(alpha=0.75, aes(x=age_in_years, y=SurfaceHoles, color=sex)) +
#   # geom_jitter(height = 0, width=0.01, alpha=0.5, aes(y=age_in_years, x=average_grade, color=sex)) +
#   scale_color_manual(values = cbbPalette[2:3], name = "Sex") +
#   labs(title="Euler Number vs. Age at Scan") + 
#   ylab("Euler Number") +
#   xlab("Age at Scan (years)") + 
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_line(color="gray"), # element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) 
# 
# print(ageEulerSex)
# c(ageHist, scannerHist, ageQcSex, eulerQcSex)
# Arrange plots
patch <- wrap_plots(demoFigs, guides = "collect", ncol=2)
png(file=fnOut,
    width=800, height=900)
print(patch) # + plot_annotation(title="Lifespan and Slip Age at Peak Volume"))
dev.off()

# LOH
# t-test sex differences: goal is not significant for demographic info
paramsToTest <- c("age_at_scan_days", "average_grade")
print("Testing sex distributions")
for (p in paramsToTest){
  testResults <- t.test(analysisDf[analysisDf$sex == "M", p], analysisDf[analysisDf$sex == "F", p])
  print(testResults$p.value)
}

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
# analysisDf <- addPrimaryScanReasonCol(analysisDf)
# reasonsTable <- summary(analysisDf$top_scan_reason_factors)
# pieLabels <- paste(names(reasonsTable), reasonsTable, sep='\n')
# pie(reasonsTable, pieLabels, col=cbbPalette, mai=c(0,0,0,0),
#     main="Top Scan Reasons for Clinical Controls\n(with sample sizes)")
