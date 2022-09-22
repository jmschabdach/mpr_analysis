library(ggplot2)
library(dplyr)
library(gamlss) #to fit model


fn <- '/Users/youngjm/Data/lifespan_growth_charts/lifespan_centile_medians.csv'
lc <- read.csv(fn)
fn <- '/Users/youngjm/Data/clip/fs6_stats/combatted_phenotypes_normalized.csv'
clipDf <- read.csv(fn)

# Prep Lifespan centiles df to match CLIP formatting for sex
lc <- lc %>%
  mutate(sex = str_replace(sex, 'Female', 'F'))
lc <- lc %>%
  mutate(sex = str_replace(sex, 'Male', 'M'))


# Add a column to clip for post conception age
clipDf$age <- clipDf$age_at_scan_days+280

# Rename the columns of the combatted clip dataframe we care about
names(clipDf)[names(clipDf) == "GMVTransformed.normalised"] <- 'GMV'
names(clipDf)[names(clipDf) == "WMVTransformed.normalised"] <- 'WMV'
names(clipDf)[names(clipDf) == "sGMVTransformed.normalised"] <- 'sGMV'
names(clipDf)[names(clipDf) == "VentriclesTransformed.normalised"] <- 'CSF'
names(clipDf)[names(clipDf) == "totalSA2Transformed.normalised"] <- 'SA'
names(clipDf)[names(clipDf) == "meanCT2Transformed.normalised"] <- 'CT'
names(clipDf)[names(clipDf) == "TCVTransformed.normalised"] <- 'TCV'

phenos <- c('GMV', 'WMV', 'sGMV', 'CSF', 'CT', 'SA', 'TCV')
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


minAgeDaysPostCon <- min(clipDf$age)
maxAgeDaysPostCon <- max(clipDf$age)
ageLimited <- unique(lc[lc$age > minAgeDaysPostCon & lc$age < maxAgeDaysPostCon, "age"])

for ( p in phenos ) {
  print(p)
  plots <- c()
  fnOut <- paste0('/Users/youngjm/Data/clip/figures/2022-09-21_clip_lifespan_correlation_', p, '.png')
  smallLc <- lc[lc$feature==p, ]
  
  # For plot 2:
  # 1. Generate GAMLSS model
  formula <- as.formula(paste(p, "pb(age_at_scan_days)", sep="~"))
  gamModelM <-gamlss(formula, 
                    family = GG, 
                    data=na.omit(clipDf[clipDf$sex == "M",]), 
                    control = gamlss.control(n.cyc = 50),
                    trace = F)
  gamModelF <-gamlss(formula, 
                     family = GG, 
                     data=na.omit(clipDf[clipDf$sex == "F",]), 
                     control = gamlss.control(n.cyc = 50),
                     trace = F)
  
  # 2. Predict the 50th centile
  clipMedM <- centiles.pred(gamModelM, type="centiles", cent=c(50), xvalues=c(sort(ageLimited)) , xname="age_at_scan_days")
  clipMedF <- centiles.pred(gamModelF, type="centiles", cent=c(50), xvalues=c(sort(ageLimited)) , xname="age_at_scan_days")
  
  # 3. Isolate Lifespan data in the correct age range for both sexes
  smallLcM <- smallLc[smallLc$sex == 'M' & smallLc$age %in% ageLimited, ]
  smallLcF <- smallLc[smallLc$sex == 'F' & smallLc$age %in% ageLimited, ]
  smallLcM <- smallLcM[order(smallLcM$age), ]
  smallLcF <- smallLcF[order(smallLcF$age), ]
  
  # 4. Calculate correlation
  rM <- cor(smallLcM$value, clipMedM$`50`)
  rF <- cor(smallLcF$value, clipMedF$`50`)
  
  # 5. Plot CLIP vs Lifespan
  # minPheno <- min(min(smallLcM$value, clipMedM$`50`),
  #                 min(smallLcF$value, clipMedF$`50`))
  # maxPheno <- max(max(smallLcM$value, clipMedM$`50`),
  #                 max(smallLcF$value, clipMedF$`50`))
  # offset <- 0.01*(maxPheno-minPheno)
  # plots[[2]] <- ggplot() +
  #   geom_point(data=clipMedF, aes(x=x, y=`50`, color=x), alpha=0.2) +
  #   geom_point(data=smallLcF, aes(x=age, y=value, color=age), alpha=0.2) +
  #   scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
  #   labs(title = paste0("Correlation Graph (Female) ", p),
  #        subtitle = paste0("(Male cor=", format(rM, nsmall=4), ", Female cor=", format(rF, nsmall=4), ")"))
   
    # plots[[2]] <- ggplot() +
    # # geom_point(aes(x=smallLcF$value, y=clipMedF$`50`, color=clipMedF$x), alpha=0.2) +
    # # geom_point(aes(x=smallLcM$value, y=clipMedM$`50`, color=clipMedM$x), alpha=0.2) +
    # geom_point(data=clipMedF, aes(x=x, y=`50`, color=x), alpha=0.2) +
    # geom_point(data=smallLcF, aes(x=age, y=value, color=age), alpha=0.2) +
    # # scale_color_manual(name = "Sex") +
    # labs(title = paste0("Correlation Graph (Female) ", p),
    #      subtitle = paste0("(Male cor=", format(rM, nsmall=4), ", Female cor=", format(rF, nsmall=4), ")"))+
    # # xlab("Lifespan 50th Centile for Phenotype") +
    # # ylab("CLIP 50th Centile for Phenotype") +
    # # xlim(minPheno-offset, maxPheno+offset) +
    # # ylim(minPheno-offset, maxPheno+offset)

  plots[[1]] <- ggplot(log='x') +
    geom_point(data=clipDf, aes(x=age, y=clipDf[,p], color=sex), alpha=0.25) +
    geom_line(data=smallLc, aes(x=age, y=value, color=sex)) +
    scale_color_manual(values = cbbPalette, name = "Sex") +
    coord_trans("log", "identity") +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(title = paste0(p))
  # # Show the plots
  # patch <- wrap_plots(plots, guides='collect')
  # png(file=fnOut,
  #     width=1200, height=600)
  # print(patch + plot_annotation(title=p))
  # dev.off()

}

