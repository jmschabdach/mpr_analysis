gc()
library(ggplot2)
library(dplyr)
library(gamlss) #to fit model
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")


fn <- '/Users/youngjm/Data/lifespan_growth_charts/lifespan_centile_medians.csv'
lc <- read.csv(fn)
fn <- '/Users/youngjm/Data/clip/fs6_stats/combatted_phenotypes_normalized.csv'
clipDf <- read.csv(fn)

# Prep Lifespan centiles df to match CLIP formatting for sex
lc <- lc %>%
  mutate(sex = str_replace(sex, 'Female', 'F'))
lc <- lc %>%
  mutate(sex = str_replace(sex, 'Male', 'M'))

# Add a column to Lifespan for log(age)
lc$logAge <- log(lc$age, 10)


# Add a column to clip for post conception age
clipDf$age <- clipDf$age_at_scan_days+280
clipDf <- addPrimaryScanReasonCol(clipDf)

# Add a column to clip for log(age)
clipDf$logAge <- log(clipDf$age, 10)

# Rename the columns of the combatted clip dataframe we care about
names(clipDf)[names(clipDf) == "GMVTransformed.normalised"] <- 'GMV'
names(clipDf)[names(clipDf) == "WMVTransformed.normalised"] <- 'WMV'
names(clipDf)[names(clipDf) == "sGMVTransformed.normalised"] <- 'sGMV'
names(clipDf)[names(clipDf) == "VentriclesTransformed.normalised"] <- 'CSF'
names(clipDf)[names(clipDf) == "totalSA2Transformed.normalised"] <- 'SA'
names(clipDf)[names(clipDf) == "meanCT2Transformed.normalised"] <- 'CT'
names(clipDf)[names(clipDf) == "TCVTransformed.normalised"] <- 'TCV'
clipDf$sex <- as.factor(clipDf$sex)
# clipDf <- clipDf[clipDf$age < max(clipDf$age), ]

phenos <- c('GMV', 'WMV', 'sGMV', 'CSF', 'CT', 'SA', 'TCV')
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

minAgeDaysPostCon <- min(clipDf$logAge)
maxAgeDaysPostCon <- max(clipDf$logAge)
ageLimited <- unique(lc[lc$logAge > minAgeDaysPostCon & lc$logAge < maxAgeDaysPostCon, "logAge"])

for ( p in phenos ) {
  print(p)
  plots <- c()
  fnOut <- paste0('/Users/youngjm/Data/clip/figures/2022-10-05_clip_lifespan_correlation_', p, '.png')
  smallLc <- lc[lc$feature==p, ]
  
  # For plot 2:
  # 1. Generate GAMLSS model
  formula <- as.formula(paste0(p, "~pb(logAge, method='GAIC', k=2) + SurfaceHoles + sex")) 
  gamModel <-gamlss(formula, 
                    family = GA, 
                    data=na.omit(clipDf), 
                    control = gamlss.control(n.cyc = 50),
                    trace = F)
  print("finished training the model")
  
  # 2. Predict phenotype values for set age range
  newData <- data.frame(logAge=sort(ageLimited),
                        SurfaceHoles=c(rep(median(clipDf$SurfaceHoles), length(ageLimited))),
                        sex=c(rep(as.factor("M"),  length(ageLimited))))
  clipPred <- predictAll(gamModel, newdata=newData)
 
  # 3. Get the smaller Lifespan data frame
  smallLc <- lc[(lc$logAge %in% ageLimited) & (lc$feature == p) & (lc$sex == 'M'), ]
  smallLc$sex <- as.factor(smallLc$sex)
  
  
  # 4. Calculate correlation
  r <- cor(smallLc$value, clipPred$mu)
   
  # 5. Plot CLIP vs Lifespan
  plots[[1]] <- ggplot(x="linear") +
    geom_point(data=clipDf, aes(x=logAge, y=clipDf[,p], color=top_scan_reason_factors), alpha=0.3) +
    geom_line(aes(x=ageLimited, y=clipPred$mu, linetype="Predicted Mu CLIP")) +
    geom_line(data=smallLc, aes(x=logAge, y=value, linetype="Lifespan 50th Centile")) +
    scale_color_manual(values = cbbPalette, name = "Reason for Scan") +
    scale_linetype_manual(values = c('solid', 'dashed'), name=" ")+
    theme(plot.title=element_text(hjust=0.5)) +
    labs(subtitle = paste0("(r=", format(r, digits=4), ")")) +
    xlab("log(Age) (days)") +
    ylab("Phenotype Value")
  
  # plots[[2]] <- ggplot(x="linear") +
  #   geom_point(data=clipDf, aes(x=age, y=clipDf[,p], color=top_scan_reason_factors), alpha=0.3) +
  #   geom_smooth(data=clipDf, aes(x=age, y=clipDf[,p])) +
  #   scale_color_manual(values = cbbPalette, name = "Reason for Scan") +
  #   theme(plot.title=element_text(hjust=0.5)) +
  #   labs(subtitle = paste0("Testing")) +
  #   xlab("Age (days)") +
  #   ylab("Phenotype Value")

  # # Show the plots
  patch <- wrap_plots(plots, guides='collect')
  png(file=fnOut,
      width=700, height=550)
  print(patch + plot_annotation(title=paste0("Comparison of Lifespan and CLIP for ", p)))
  dev.off()

}

