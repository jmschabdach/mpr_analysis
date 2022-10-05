gc()
source("/Users/youngjm/Projects/mpr_analysis/r/lib_mpr_analysis.r")

## Step 7: Prepare regional phenotype data -----------------------

# Load the data
fn <- '/Users/youngjm/Data/clip/fs6_stats/original_phenotypes_singleScanPerSubject.csv'
t <- "Age at Peak of Original Regional Phenotypes"
fnOut <- '/Users/youngjm/Data/clip/figures/2022-09-19_hemisphere_original_age_at_peak.png'

brainDf <- read.csv(fn)
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
  # script currently barfing here, subscript out of bounds
  print(parsedRegionalPhenotypes[[i]])
  print(grayVolLhCols[[i]])
  print(grayVolRhCols[[i]])
  brainDf[[parsedRegionalPhenotypes[[i]]]] <- (brainDf[[grayVolLhCols[[i]]]] + brainDf[[grayVolRhCols[[i]]]])/2
}

globalPhenoCols <- c('TotalGrayVol', 'CerebralWhiteMatterVol', 'SubCortGrayVol',
                     'eTIV', 'VentricleVolume', 'CorticalSurfaceArea', 
                     'MeanCorticalThickness', 'TCV')
regionalMetaCols <- setdiff(brainDfCols, c(globalPhenoCols, grayVolCols))
greyVolDf <- brainDf[append(regionalMetaCols, parsedRegionalPhenotypes)]

## Results 4: Make regional phenotype plots for age at peak --------------------
##
# Make the composite plot and tables for a given dataframe/set of phenotypes
# @param df A dataframe where rows are scans and columns are phenotypes and features
# @param phenotypes A vector of phenotypes to analyze
# @param colNames A vector of cleaned names to use for the table columns
# @param rowNames A vector of cleaned names to use for the table rows
# @param title A string to differentiate the figures for this dataframe vs. other dataframes
generateHemispherePlots <- function(df, phenotypes, phenotypesNoWS, title, sex=""){
  # Initialize variables
  # Set up empty variables to generate table later
  ageAtPeak <- c()
  
  # for phenotype in phenotypes...
  for (phenotype in phenotypesNoWS){
    print(phenotype)
    if (sex == "M") {
      print("M")
      modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
                               sex="M",
                               top_scan_reason_factors='headaches')
      gamm <- createGammNoSex(df[df$sex == "M", ], phenotype)
      gammPreds <- predict_gam(gamm, values = modelFixedValues)
    } 
    else if (sex == "F"){
      print("F")
      modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
                               sex="F",
                               top_scan_reason_factors='headaches')
      gamm <- createGammNoSex(df[df$sex == "F", ], phenotype)
      gammPreds <- predict_gam(gamm, values = modelFixedValues)
    } 
    else {
      modelFixedValues <- list(SurfaceHoles = mean(df$SurfaceHoles), 
                               sex='M',
                               top_scan_reason_factors='headaches')
      # Generate GAMs using the formula that incorporates scan reason
      gamm <- createMainGamm(df, phenotype)
      print("b")
      # Predict on the GAMs for the actual data
      gammPreds <- predict_gam(gamm, values = modelFixedValues)
      print("c")
    }
    
    # Get the age at peak
    ageAtPeak <- append(ageAtPeak, c(getAgeAtPeak(gammPreds)))
  }
  
  # Convert the lists into a dataframe
  results = as.data.frame(cbind(region=phenotypes, 
                                feat=phenotypesNoWS,
                                em=as.numeric(ageAtPeak)),
                          stringsAsFactors=F)
  
  print(results)
  
  # Plot the dataframe using ggseg
  p <- results %>% 
    ggseg(mapping=aes(fill=as.numeric(em)),
          hemisphere='left') +
    labs(title = "Age at Peak (years)", legend="Age") +
    theme(axis.title = element_blank()) +
    scale_fill_gradient(low = "blue", high = "red", na.value = NA) 
  
  grid.arrange(p)
  return(results)
}

regionalPeaks <- generateHemispherePlots(greyVolDf, regionalPhenotypes, parsedRegionalPhenotypes, '')
regionalPeaksM <- generateHemispherePlots(greyVolDf, regionalPhenotypes, parsedRegionalPhenotypes, '', "M")
regionalPeaksF <- generateHemispherePlots(greyVolDf, regionalPhenotypes, parsedRegionalPhenotypes, '', "F")
