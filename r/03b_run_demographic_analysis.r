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

# Filepaths to change
fn <- '/Users/youngjm/Data/slip/fs6_stats/original_phenotypes_demographics.csv'
fnOut <- '/Users/youngjm/Data/slip/figures/2022-11-07_demographics_figure.tiff'

analysisDf <- read.csv(fn)

# Convert some columns to factors
toFactor <- c('sex', 'fs_version', 'MagneticFieldStrength', 'scanner_id', 'scan_reason_primary')
analysisDf[toFactor] <- lapply(analysisDf[toFactor], factor)

demoFigs <- c()

analysisDf$scanner_id <-gsub("Scanner ","",as.character(analysisDf$scanner_id))

# Make a bar chart for age at scan
demoFigs[[1]] <- ggplot(analysisDf, aes(x=age_in_years, fill=sex)) +
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
results <- cor.test(analysisDf$average_grade, analysisDf$age_in_years, method="p", conf.level=0.95)
r <- results$estimate
pval <- results$p.value
ciLow <- results$conf.int[1]
ciHigh <- results$conf.int[2]
demoFigs[[3]] <- ggplot(analysisDf) +
  geom_hline(yintercept=1.0, color="black", alpha=0.5) +
  # geom_point(alpha=0.75, aes(x=age_in_years, y=average_grade, color=sex)) +
  geom_jitter(height = 0.025, width=0.01, alpha=0.5, aes(x=age_in_years, y=average_grade, color=sex)) +
  scale_color_manual(values = cbbPalette[2:3], name = "Sex") +
  labs(title="Average QC Rating vs. Age at Scan",
       subtitle=paste0("(r=", format(r, digits=3), ", p=", format(pval, digits=3), 
                       " (95% CI: ", format(ciLow, digits=3),", ",format(ciHigh, digits=3),  "))")) +
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
results <- cor.test(analysisDf$SurfaceHoles, analysisDf$average_grade, method="p", conf.level=0.95)
r <- results$estimate
pval <- results$p.value
ciLow <- results$conf.int[1]
ciHigh <- results$conf.int[2]
demoFigs[[5]] <- ggplot(analysisDf) +
  geom_hline(yintercept=1.0, color="black", alpha=0.5) +
  geom_jitter(height = 0.1, width=0.1, alpha=0.5, aes(x=SurfaceHoles, y=average_grade, color=sex)) +
  scale_color_manual(values = cbbPalette[2:3], name = "Sex") +
  labs(title=paste0("Average QC Rating vs. Euler Number"),
       subtitle=paste0("(r=", format(r, digits=3), ", p=", format(pval, digits=3), 
                       " (95% CI: ", format(ciLow, digits=3),", ",format(ciHigh, digits=3),  "))")) + 
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

# Arrange plots
patch <- wrap_plots(demoFigs, guides = "collect", ncol=2)
tiff(file=fnOut,
    width=6000, height=6000, res=600)
print(patch) 
dev.off()

patch <- wrap_plots(demoFigs[[3]] + demoFigs[[4]] + demoFigs[[5]] + demoFigs[[6]], 
                    guides="collect")
tiff(file="/Users/youngjm/Documents/80 Communications/83 Papers/83.02 CLIP Methods/Radiology Submission Figures/Supplemental Figure 01 additional demographics.tiff",
     width=6000, height=4000, res=600)
print(patch)
dev.off()

# LOH
# t-test sex differences: goal is not significant for demographic info
paramsToTest <- c("age_at_scan_days", "average_grade")
print("Testing sex distributions")
for (p in paramsToTest){
  testResults <- t.test(analysisDf[analysisDf$sex == "M", p], analysisDf[analysisDf$sex == "F", p])
  print(p)
  print(testResults$p.value)
  print(testResults$conf.int)
  print(testResults$method)
  print(" ")
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

