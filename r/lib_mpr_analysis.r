# #-------------------------------------------------------------------------------
# # FUNCTION DEFINITIONS
# #-------------------------------------------------------------------------------
# 
# ## Demographic plot functions:
# generateScannerDistributionPlot <- function(dataDf, sex, colorPalette){
#   # Make a set of labels for each scanner
#   scannerLabels <- str_pad(c(1:length(levels(dataDf$scanner_id))), 2, pad = "0")
#   scannerLabels <- paste('Scanner', scannerLabels, sep='\n')
#   print(scannerLabels)
#   
#   # Set the title string
#   if (sex == 'M'){
#     title = "Scanner Distribution (Male)"
#     dataDf <- dataDf[dataDf$sex == 'M', ]
#   } else {
#     title = "Scanner Distribution (Female)"
#     dataDf <- dataDf[dataDf$sex == 'F', ]
#   }
#   
#   # dataDf <- data.frame(table(dataDf))
# 
#   newPlot <- ggplot(data=dataDf, aes(x=scanner_id, fill=as.factor(rawdata_image_grade))) +
#     geom_bar(position="stack") +
#     scale_fill_manual(values = qcColors, drop=FALSE) +
#     scale_x_discrete(labels=scannerLabels, drop=FALSE) +
#     labs(title = title,
#          x = "Scanner ID",
#          y = "# Scans per Scanner",
#          fill = "Image Quality") +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   return(newPlot)
# }
# 
# generateAgeDistributionPlot <- function(dataDf, sex, colorPalette){
#   # Get the max age
#   xmax <- max(dataDf$age_in_years)+1
#   print(xmax)
#   
#   # Set the title string
#   if (sex == 'M'){
#     title = "Age at Scan (Male)"
#     dataDf <- dataDf[dataDf$sex == 'M', ]
#   } else {
#     title = "Age at Scan (Female)"
#     dataDf <- dataDf[dataDf$sex == 'F', ]
#   }
#   
#   newPlot <- ggplot(data=dataDf, aes(x=age_in_years, fill=as.factor(rawdata_image_grade))) +
#     geom_histogram(position=position_stack(), binwidth = 0.5) +
#     scale_fill_manual(values = qcColors, drop=FALSE) + 
#     xlim(0, xmax) +
#     labs(title = title,
#          x = 'Patient Age (Years)',
#          y = '# Patients',
#          fill="Image Quality") +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   return(newPlot)
# }
# 
# generateEulerQcDistributionPlot <- function(dataDf, sex, colorPalette){
#   # Set the title string
#   if (sex == 'M'){
#     title = "Euler Number vs. QC Rating (Male)"
#     dataDf <- dataDf[dataDf$sex == 'M', ]
#   } else {
#     title = "Euler Number vs. QC Rating (Female)"
#     dataDf <- dataDf[dataDf$sex == 'F', ]
#   }
#   
#   # Make plots for QC
#   pQcSurfaceHoles <- dataDf %>%
#     ggplot(aes(x=as.factor(rawdata_image_grade), y=SurfaceHoles)) +
#     geom_violin(aes(fill=as.factor(rawdata_image_grade))) +
#     geom_boxplot(width=0.1)+ 
#     scale_fill_manual(values=qcColors) +
#     labs(title = title,
#          x = 'Image QC Rating Group',
#          y = 'Number of Surface Holes (count)',
#          fill="Image Quality") +
#     theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
#   
#   # Estimate the effect?
# }
# 
# generateAgeQcDistributionPlot <- function(dataDf, sex, colorPalette){
#   
#   # Set the title string
#   if (sex == 'M'){
#     title = "Age vs. QC Rating (Male)"
#     dataDf <- dataDf[dataDf$sex == 'M', ]
#   } else {
#     title = "Age vs. QC Rating (Female)"
#     dataDf <- dataDf[dataDf$sex == 'F', ]
#   }
#   
#   pQcAge <- dataDf %>%
#     ggplot(aes(x=as.factor(rawdata_image_grade), y=age_in_years)) +
#     geom_violin(aes(fill=as.factor(rawdata_image_grade))) +
#     geom_boxplot(width=0.1)+ 
#     scale_fill_manual(values=qcColors)+
#     labs(title = title,
#          x = 'Image QC Rating Group',
#          y = 'Age at Scan (Days)',
#          fill="Image Quality") +
#     theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
# }
# 
# ##
# # Given a dataframe and a column header (imaging phenotype), create a GAMM
# # (GAMM = Generalized Additive Mixed Model)
# # @param df A dataframe where the rows are scans and columns contain phenotypes
# # @param measure Name of a column in the dataframe
# # @return mixedModel The mixed model for the given formula and specificed measure
# createMainGamm <- function(df, measure) {
#   # also try k=10
#   formula <- as.formula(paste(measure, "s(log(age_in_years), k=3) + 
#                       SurfaceHoles +
#                       sex", sep="~"))
#   
#   # Make a basic linear model accounting for age and surface holes for the given data frame
#   mixedModel <- gam(formula,
#                     data = df,
#                     gamma=1)
#   # Return the model
#   return(mixedModel)
# }
# 
# createGammLinearAge <- function(df, measure) {
#   formula <- as.formula(paste(measure, "log(age_in_years) +
#                       SurfaceHoles +
#                       # top_scan_reason_factors +
#                       sex", sep="~"))
#   
#   # modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
#   #                          sex='M', 
#   #                          top_scan_reason_factors='headaches')
#   
#   # Make a basic linear model accounting for age and surface holes for the given data frame
#   mixedModel <- gamm(formula,
#                      random = list(scanner_id=~1),
#                      data = df)
#   # Return the model
#   return(mixedModel)
# }
# 
# createGammScanReasonsAge <- function(df, measure) {
#   formula <- as.formula(paste(measure, "s(log(age_in_years), fx=T) +
#                       SurfaceHoles +
#                       top_scan_reason_factors +
#                       sex", sep="~"))
#   
#   # modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
#   #                          sex='M', 
#   #                          top_scan_reason_factors='headaches')
#   
#   # Make a basic linear model accounting for age and surface holes for the given data frame
#   mixedModel <- gamm(formula,
#                      random = list(scanner_id=~1),
#                      data = df)
#   # Return the model
#   return(mixedModel)
# }
# 
# createGammNoSex <- function(df, measure) {
#   formula <- as.formula(paste(measure, "s(log(age_in_years), k=10) +
#                       SurfaceHoles", sep="~"))
#   
#   # Make a basic linear model accounting for age and surface holes for the given data frame
#   mixedModel <- gam(formula,
#                     data = df,
#                     gamma=1)
#   # Return the model
#   return(mixedModel)
# }
# 
# ##
# # Make a scatter plot of a specific measure for the original data. Include the
# # confidence intervals that come from the output of the GAMM predictions
# # @param predictions Output from predict_gam() for the measure and its gam
# # @param origData A dataframe containing the original data
# # @param measure A column of the original dataframe
# # @param measureTitle A string used to describe the measure in the plot
# # @returns plot01 A ggplot object with scatter points and a confidence interval
# generatePlotScatterWithCI <- function(predictions, origData, measure, measureTitle){
#   # Colorblind palette with black:
#   cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#   
#   peakAge <- getAgeAtPeak(predictions)
#   
#   plot01 <- predictions %>%
#     ggplot(aes(x=age_in_years, y=fit), log='x') +
#     geom_point(data=origData, alpha=0.5, aes(x=age_in_years, y=measure, color=sex)) +
#     geom_smooth_ci() +
#     geom_vline(xintercept=peakAge, linetype='dashed', color='#c23400') +
#     scale_color_manual(values = cbbPalette, name = "Sex") +
#     # scale_x_continuous(trans='log10') +
#     theme(axis.title = element_blank(), plot.title=element_text(hjust=0.5)) +
#     labs(title = measureTitle)
#   
#   print(plot01)
#   
#   return(plot01)
# }
# 
# generateDiagnosisPlotCI <- function(predictions, origData, measure, measureTitle){
#   
#   dxColors <- c("#000000", "#E69F00")
#   
#   plot02 <- predictions %>%
#     ggplot(aes(x=age_in_years, y=fit)) +
#     geom_smooth_ci() +
#     scale_color_manual(values = dxColors) +
#     # scale_x_continuous(trans='log10') +
#     theme(axis.title = element_blank()) 
#   # labs(title = "Growth Trajectory Over Lifespan") #paste(measureTitle, ": Growth Trajectories"))
#   
#   return(plot02)
# }
# 
# 
# ##
# # Make a table showing the statistical significance of nonlinear parameters on image phenotype
# # @param aGam A gam object
# # @param nTests A number indicating how many tests are being performed/should be corrected for
# # @returns paraT A table object to be displayed using gridExtra
# generateScanReasonsParametricTable <- function(aGam){
#   subTab <- signif(summary(aGam)$s.table, 3)
#   rownames(subTab) <- c('Age')
#   
#   # Add the p.adjust
#   # subTab <- cbind(subTab, p.adjust(subTab[,"p-value"], n = nTests))
#   colnames(subTab) <- c('edf', 'Ref.df', 'F', 'p-value')
#   cols <- c("p-value")
#   subTab <- subTab[,cols]
#   
#   return(subTab)
# }
# 
# 
# 
# ##
# # Make a table showing the statistical significance of linear parameters on image phenotype
# # @param aGam A gam object
# # @param nTests A number indicating how many tests are being performed/should be corrected for
# # @returns linT A table object to be displayed using gridExtra
# generateScanReasonsLinearTable <- function(aGam){
#   subTab <- signif(summary(aGam)$p.table, 3)
#   
#   if ("sexM" %in% dimnames(subTab)[[1]] | "sexF" %in% dimnames(subTab)[[1]]){
#     newRowNames <- c('Surface Holes',
#                      'Sex')
#   }
#   else {
#     newRowNames <- c('SurfaceHoles')
#   }
#   
#   rows <- rownames(subTab)[2:length(rownames(subTab))]
#   subTab <- subset(subTab, rownames(subTab) %in% rows)
#   rownames(subTab) <- newRowNames
#   
#   # Rename and trim the columns
#   colnames(subTab) <- c('beta', 'std err', 't-value', 'p-value')
#   cols <- c('p-value')
#   subTab <- subTab[, colnames(subTab) %in% cols]
#   
#   if (length(newRowNames) == 1){
#     subTab <- data.frame(subTab, row.names = cols)
#     colnames(subTab) <- newRowNames
#   }
#   
#   return(subTab)
# }
# 
# ##
# # Make the composite plot and tables for a given dataframe/set of phenotypes
# # @param df A dataframe where rows are scans and columns are phenotypes and features
# # @param phenotypes A vector of phenotypes to analyze
# # @param colNames A vector of cleaned names to use for the table columns
# # @param rowNames A vector of cleaned names to use for the table rows
# # @param title A string to differentiate the figures for this dataframe vs. other dataframes
# generatePlotsAndTablesForDataset <- function(df, phenotypes, colNames, rowNames, title, sex=""){
#   # Initialize variables
#   # Set up empty variables to generate table later
#   ageAtPeak <- list()
#   predictionPlots <- list()
#   tableValues <- list()
#   tableBetas <- list()
#   
#   # for phenotype in phenotypes...
#   for (i in 1:length(phenotypes)){
#     phenotype <- phenotypes[[i]]
#     col <- colNames[i]
#     if (sex == "M"){
#       modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
#                                top_scan_reason_factors='headaches')
#       gamm <- createGammNoSex(df[df$sex == "M", ], phenotype)
#       gammPreds <- predict_gam(gamm, values = modelFixedValues)
#     } 
#     else if (sex == "F") {
#       modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
#                                top_scan_reason_factors='headaches')
#       gamm <- createGammNoSex(df[df$sex == "F", ], phenotype)
#       gammPreds <- predict_gam(gamm, values = modelFixedValues)
#     }
#     else {
#       modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
#                                sex='M',
#                                top_scan_reason_factors='headaches')
#       # Generate GAMs using the formula that incorporates scan reason
#       gamm <- createMainGamm(df, phenotype)
#       # Predict on the GAMs for the actual data
#       gammPreds <- predict_gam(gamm, values = modelFixedValues)
#     }
#     
#     # Get the age at peak
#     ageAtPeak[[phenotype]] <- getAgeAtPeak(gammPreds)
#     # Generate scatter plots with confidence intervals
#     plotScatterScanReason <- generatePlotScatterWithCI(gammPreds, df, df[ , phenotype], col)
#     predictionPlots[[phenotype]] <- plotScatterScanReason
#     
#     # tParam <- generateScanReasonsParametricTable(gammScanReason$gam)
#     # tLin <- generateScanReasonsLinearTable(gammScanReason$gam)
#     print("a")
#     tParam <- generateScanReasonsParametricTable(gamm)
#     print("b")
#     tLin <- generateScanReasonsLinearTable(gamm)
#     tableValues[[phenotype]] <- append(tLin, list(Age = tParam))
#     
#     # Get normalized betas
#     normDf <- data.frame(df)
#     normDf[[phenotype]] <- scale(normDf[[phenotype]])
#     normDf$SurfaceHoles <- scale(normDf$SurfaceHoles)
#     normDf$age_in_years <- scale(normDf$age_in_years)
#     gammNorm <- createMainGamm(normDf, phenotype)
#     # tableBetas[[phenotype]] <- getBetasFromGamSummary(gammNorm$gam)
#     tableBetas[[phenotype]] <- getBetasFromGamSummary(gammNorm)
#   }
#   
#   # Plot all scatter/CI plots in 1 figure
#   patchwork <- wrap_plots(predictionPlots, nrow=1, guides="collect")
#   print(patchwork + plot_annotation(title=paste("Trajectories of", title)))
#   # Make table: p-values of factor/phenotype
#   firstStep <- lapply(tableValues, unlist)
#   pvals <- as.data.frame(firstStep, stringsAsFactors = F)
#   colnames(pvals) <- colNames
#   rownames(pvals) <- rowNames
#   # t1 <- gridExtra::tableGrob(secondStep)
#   # gridExtra::grid.arrange(top=paste("p-values of", title), t1) 
#   print("woo 1")
#   
#   # Make table: betas
#   firstStep <- lapply(tableBetas, unlist)
#   betaVals <- as.data.frame(firstStep, stringsAsFactors = F)
#   colnames(betaVals) <- tableColumns
#   rownames(betaVals) <- tableRows[1:length(tableRows)-1]
#   # t2 <- gridExtra::tableGrob(secondStep)
#   # gridExtra::grid.arrange(top=paste("Betas of", title), t1)
#   print("woo 2")
#   
#   output <- list()
#   output$predictionPlots <- predictionPlots
#   output$pvalsTable <- pvals
#   output$betasTable <- betaVals
#   output$ageAtPeak <- ageAtPeak
#   
#   return(output)
# }
# 
# ##
# # Make the composite plot and tables for a given dataframe/set of phenotypes
# # @param df A dataframe where rows are scans and columns are phenotypes and features
# # @param phenotypes A vector of phenotypes to analyze
# # @param colNames A vector of cleaned names to use for the table columns
# # @param rowNames A vector of cleaned names to use for the table rows
# # @param title A string to differentiate the figures for this dataframe vs. other dataframes
# generateCentilePlots <- function(df, phenotypes, cols, title){
#   # Initialize variables
#   centilePlots <- list()
#   cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#   
#   # for phenotype in phenotypes...
#   for (i in 1:length(phenotypes)){
#     phenotype <- phenotypes[[i]]
#     # Generate scatter plots with confidence intervals
#     plot01 <- df %>%
#       ggplot(aes(x=age_in_years, y=.data[[phenotype]]), log='x') +
#       geom_point(data=df, alpha=0.5, aes(color=sex)) +
#       scale_color_manual(values = cbbPalette, name = "Sex") +
#       theme(axis.title = element_blank(), plot.title=element_text(hjust=0.5)) +
#       labs(title = cols[[i]])
#     
#     centilePlots[[phenotype]] <- plot01
#   }
#   
#   # Plot all scatter/CI plots in 1 figure
#   patchwork <- wrap_plots(centilePlots, nrow=1, guides="collect")
#   print(patchwork + plot_annotation(title=paste("Centiles of", title)))
#   
#   return(centilePlots)
# }
# 
# 
# generatePlotsAndTablesLinearAge <- function(df, phenotypes, colNames, rowNames, title){
#   # Initialize variables
#   # Set up empty variables to generate table later
#   ageAtPeak <- list()
#   predictionPlots <- list()
#   tableValues <- list()
#   tableBetas <- list()
#   
#   modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
#                            sex='M', 
#                            top_scan_reason_factors='headaches')
#   
#   # for phenotype in phenotypes...
#   for (phenotype in phenotypes){
#     # Generate GAMs using the formula that incorporates scan reason
#     gammScanReason <- createGammLinearAge(df, phenotype)
#     
#     # Predict on the GAMs for the actual data
#     gammScanReasonPreds <- predict_gam(gammScanReason$gam, values = modelFixedValues)
#     
#     # Get the age at peak
#     ageAtPeak[[phenotype]] <- getAgeAtPeak(gammScanReasonPreds)
#     
#     # Generate scatter plots with confidence intervals
#     plotScatterScanReason <- generatePlotScatterWithCI(gammScanReasonPreds, df, df[ , phenotype], phenotype)
#     predictionPlots[[phenotype]] <- plotScatterScanReason
#     
#     tLin <- signif(summary(gammScanReason$gam)$p.table, 3)
#     colnames(tLin) <- c('beta', 'std err', 't-value', 'p-value')
#     tableValues[[phenotype]] <- tLin
#   }
#   
#   # Plot all scatter/CI plots in 1 figure
#   patchwork <- wrap_plots(predictionPlots, nrow=1, guides="collect")
#   print(patchwork + plot_annotation(title=paste("Trajectories of", title)))
#   
#   return(tableValues)
# }
# 
# ##
# # Given a set of data, a phenotype, and an age range to predict over, 
# #   1. Generate a GAMLSS model in the GG family with nu set to the identity function
# #   2. Predict the phenotype values for the specified age range
# # @param y String specifying the  phenotype to predict
# # @param df A dataframe containing columns for the phenotype, the log of the postconception age (logAge), the Euler number (SurfaceHoles), and sex as a factor
# # @param ageForPred A list of age points to perform predictions at
# # @return predictedModel The trained GAMLSS model predictions on the age range
# predictGamlssModel <- function(y, df, ageForPred) {
#   # 1. Generate GAMLSS models
#   formula <- as.formula(paste0(y, "~fp(logAge, npoly=3) + SurfaceHoles + sex")) 
#   gamModel <-gamlss(formula = formula, 
#                     sigma.formula = formula,
#                     nu.formula = as.formula(paste0(y, "~1")),
#                     family = GG, 
#                     data=na.omit(df), 
#                     control = gamlss.control(n.cyc = 200),  # lifespan
#                     trace = F)
#   print("Model trained.")
#   
#   # 2. Predict phenotype values for set age range
#   newData <- data.frame(logAge=sort(ageForPred),  # possibly put a thing in there to calculate the limited age range for me
#                         SurfaceHoles=c(rep(median(df$SurfaceHoles), length(ageForPred))),
#                         sex=c(rep(as.factor("M"),  length(ageForPred))))
#   predictedModel <- predictAll(gamModel, newdata=newData)
#   print("Predictions generated for specified age range.")
#   
#   return(predictedModel)
# }

##
#
#
#
buildFreeSurferGAMM <- function(p, df) {
  # 1. Generate GAMLSS models
  formula <- as.formula(paste0(p, "~fp(logAge, npoly=3) + SurfaceHoles + sex - 1"))
  gamModel <-gamlss(formula = formula,
                    sigma.formula = formula,
                    nu.formula = as.formula(paste0(p, "~1")),
                    family = GG,
                    data = na.omit(df),
                    control = gamlss.control(n.cyc = 200),  # Per lifespan nature paper
                    trace = F)
  print("finished training the models")
  return(gamModel)
}

predictCentilesForAgeRange <- function(gamModel, ageRange, euler=0, cent=0.5){
  if (euler == 0){
    newDataM <- data.frame(logAge=sort(ageRange),
                           sex=c(rep(as.factor("M"),  length(ageRange))))
    newDataF <- data.frame(logAge=sort(ageRange),
                           sex=c(rep(as.factor("F"),  length(ageRange))))
  } else  {
    print("euler!")
    newDataM <- data.frame(logAge=sort(ageRange),
                           SurfaceHoles=c(rep(euler, length(ageRange))),
                           sex=c(rep(as.factor("M"),  length(ageRange))))
    
    newDataF <- data.frame(logAge=sort(ageRange),
                           SurfaceHoles=c(rep(euler, length(ageRange))),
                           sex=c(rep(as.factor("F"),  length(ageRange))))
  } 
  
  print(cent)
  
  # Predict phenotype values for set age range for each sex
  gammModelM <- predictAll(gamModel, newdata=newDataM)
  gammModelF <- predictAll(gamModel, newdata=newDataF)
  
  # Calculate the `cent`th centiles for the sex models
  phenoMedianPredsM <- qGG(c(cent), 
                           mu=gammModelM$mu, 
                           sigma=gammModelM$sigma, 
                           nu=gammModelM$nu)
  
  phenoMedianPredsF <- qGG(c(cent), 
                           mu=gammModelF$mu, 
                           sigma=gammModelF$sigma, 
                           nu=gammModelF$nu)
  # Average the two calculated centile curves
  phenoMedianPreds <- (phenoMedianPredsF + phenoMedianPredsM)/2
  
  # Return the calculated centile curve
  return(phenoMedianPreds)
}


##
# Add a column to the dataframe - CONFIRMED
addPrimaryScanReasonCol <- function(df){
  # if (!'top_scan_reason_factors' %in% colnames(df)){
  # get the top 6 primary reasons
  topScanReasons <- names(sort(table(df$scan_reason_primary), decreasing=TRUE)[1:4])
  print(topScanReasons)
    
  # build a column with these top 5 primary reasons and a generous other category
  df <- mutate(df, top_scan_reason_factors = if_else(is.element(scan_reason_primary, topScanReasons),
                                                     paste(scan_reason_primary),
                                                     "other"))
  # }
  
  df$top_scan_reason_factors <- as.factor(df$top_scan_reason_factors)
  # Put the "other" category first
  df$top_scan_reason_factors <- relevel(df$top_scan_reason_factors, "other")
  
  return(df)
}

# PULLED FROM SCRIPTS
# Calculating the centile for a subject
calculatePhenotypeCentile <- function(model, measuredPhenotypeValue, logAge, surfaceHoles, sex){
  centileDistribution <- 1:9999/10000
  centiles <- c()
  for (i in 1:length(measuredPhenotypeValue)){
    newData <- data.frame(logAge=logAge[[i]],
                          SurfaceHoles=surfaceHoles[[i]],
                          sex=sex[[i]])
    predModel <- predictAll(model, newdata=newData)
    expectedPhenotypeValue <- qGG(centileDistribution, mu=predModel$mu, sigma=predModel$sigma, nu=predModel$nu)
    centiles[i] <- centileDistribution[which.min(abs(measuredPhenotypeValue[[i]] - expectedPhenotypeValue))]
  }
  return(centiles)
}

calculatePhenotypeCentileSynthSeg <- function(model, measuredPhenotypeValue, logAge, surfaceHoles, sex){
  centileDistribution <- 1:9999/10000
  centiles <- c()
  for (i in 1:length(measuredPhenotypeValue)){
    newData <- data.frame(logAge=logAge[[i]],
                          sex=sex[[i]])
    predModel <- predictAll(model, newdata=newData)
    expectedPhenotypeValue <- qGG(centileDistribution, mu=predModel$mu, sigma=predModel$sigma, nu=predModel$nu)
    centiles[i] <- centileDistribution[which.min(abs(measuredPhenotypeValue[[i]] - expectedPhenotypeValue))]
  }
  return(centiles)
}

convertAgeToYears <- function(ages){
  for (i in names(ages)){
    print(i)
    print((ages[[i]]-280)/365.25)
  }
}

calculatePipelineQc <- function(metrics1, metrics2){
  qc <- c()
  for (i in c(1:length(metrics1))) {
    qc <- append(qc, abs(metrics1[i]-metrics2[i])/((metrics1[i]+metrics2[i])/2))
  }
  return(qc)
}

prepForCombat <- function(df, fnBase){
  # Move metadata to the front of the dataframe
  cols <- colnames(df)
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",         
                "proc_ord_year",
                "MagneticFieldStrength", "scanner_id", #"confirm_neurofibromatosis", 
                "rawdata_image_grade", 'average_grade', 'fs_version', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories", 'SurfaceHoles')
  phenoCols <- setdiff(cols, metaCols)
  globalPhenoCols <- c('GMV', 'WMV', 'sGMV', 'CSF', 'TCV')
  nonGlobalCols <- setdiff(phenoCols, globalPhenoCols)
  
  regionalPhenoCols <- c()
  for (c in nonGlobalCols){
    if ((grepl('lh_', c) | grepl('rh_', c)) & grepl('_grayVol', c)){
      regionalPhenoCols <- append(regionalPhenoCols, c)
    }
  }
  
  newCols <- c(metaCols, globalPhenoCols, regionalPhenoCols)
  # print(setdiff(newCols, colnames(df)))
  df <- df[, newCols]
  
  # Drop any scans with NAs
  df <- df[complete.cases(df), ]
  print(colnames(df))
  
  # Identify the scan_id + phenotypes to combat
  toCombat <- df[, c('scan_id', globalPhenoCols, regionalPhenoCols)]
  
  # Identify covariates to harmonize on/protect
  covars <- data.frame(SITE = as.character(df$scanner_id), # <-- harmonize on the first column
                       log_age_in_years = log(df$age_in_years), # <-- protect all following columns
                       sex = df$sex,
                       SurfaceHoles = df$SurfaceHoles,
                       reason = df$top_scan_reason_factors)
  
  # Save the data and covars dataframes to csvs
  write.csv(toCombat,paste0(fnBase, "_phenotypes.csv"), row.names = FALSE)
  write.csv(covars,paste0(fnBase, "_covariates.csv"), row.names = FALSE)
}

prepForCombatSynthSeg <- function(df, fnBase){
  # Move metadata to the front of the dataframe
  cols <- colnames(df)
  print(cols)
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",     
                "proc_ord_year",
                "MagneticFieldStrength", "scanner_id", #"confirm_neurofibromatosis", 
                "rawdata_image_grade", 'average_grade', 'fs_version', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories")
  phenoCols <- setdiff(cols, metaCols)
  globalPhenoCols <- c('GMV', 'WMV', 'sGMV', 'CSF', 'TCV')
  nonGlobalCols <- setdiff(phenoCols, globalPhenoCols)
  
  regionalPhenoCols <- c()
  
  newCols <- c(metaCols, globalPhenoCols, regionalPhenoCols)
  df <- df[, newCols]
  
  # Drop any scans with NAs
  df <- df[complete.cases(df), ]
  
  # Identify the scan_id + phenotypes to combat
  toCombat <- df[, c('scan_id', globalPhenoCols, regionalPhenoCols)]
  
  # Identify covariates to harmonize on/protect
  covars <- data.frame(SITE = as.character(df$scanner_id), # <-- harmonize on the first column
                       log_age_in_years = log(df$age_in_years), # <-- protect all following columns
                       sex = df$sex,
                       # SurfaceHoles = df$SurfaceHoles,
                       reason = df$top_scan_reason_factors)
  
  # Save the data and covars dataframes to csvs
  write.csv(toCombat,paste0(fnBase, "_phenotypes.csv"), row.names = FALSE)
  write.csv(covars,paste0(fnBase, "_covariates.csv"), row.names = FALSE)
}

loadCombattedData <- function(df, fn){
  metaCols <- c('patient_id', "scan_id", "age_at_scan_days", 'age_in_years', "sex",      
                "proc_ord_year", 
                "MagneticFieldStrength", "scanner_id", #"confirm_neurofibromatosis", 
                "rawdata_image_grade", 'fs_version', 'top_scan_reason_factors', 
                "scan_reason_primary", "scan_reason_categories")
  if ("SurfaceHoles" %in% colnames(df)){
    metaCols <- append(metaCols, "SurfaceHoles")
  } 
  combattedDf <- read.csv(fn)
  combattedDf$scan_id <- as.factor(combattedDf$scan_id)
  df$scan_id <- as.factor(df$scan_id)
  print(dim(df))
  print(dim(combattedDf))
  # Drop rows from the dataframe if they are not in the combattedDf
  df <- df[ df$scan_id %in% combattedDf$scan_id, ]
  print(dim(df))
  combattedDf <- combattedDf[ combattedDf$scan_id %in% df$scan_id, ]
  print(dim(combattedDf))
  # Drop scan_ids column from the combattedDf (duplicate)
  combattedDf <- merge(combattedDf, df[, metaCols], by='scan_id')
  print(dim(combattedDf))
  # Add the metadata back in
  return(combattedDf)
}


