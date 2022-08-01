from neuroHarmonize import harmonizationLearn
import pandas as pd
import numpy as np

# Load the data
phenoFn = '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_phenotypes_toCombat.csv'
covFn = '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_covariates_toCombat.csv'

phenoDf = pd.read_csv(phenoFn)
covDf = pd.read_csv(covFn)

# REMOVE THIS FOR OTHER DATA, ISSUE WITH A SINGLE DATASET 
phenoDf = phenoDf.drop(index=[379])
covDf = covDf.drop(index=[379])

# Pull out the scan ids from the dataframe
scanIds = phenoDf['scan_id']
data = np.array(phenoDf.drop(columns=['scan_id']))

# Specify categorical columns
catCols = ['sex', 'reason']

# Set up the covariates
covDf = pd.get_dummies(covDf, columns=catCols, drop_first=True)

# Run ComBat
# -  Using the return_s_data flag returns a third arg: combatted data with covariate effects not preserved
model, data_combatted_covpreserved, data_combatted_nocov= harmonizationLearn(data, covDf, return_s_data=True)

# Convert the data (numpy.array) into dataframes (pd.DataFrame)
dfCombattedCovPreserved = pd.DataFrame(columns=list(phenoDf.drop(columns=['scan_id'])), data=data_combatted_covpreserved)
dfCombattedCovRemoved = pd.DataFrame(columns=list(phenoDf.drop(columns=['scan_id'])), data=data_combatted_nocov)

# Add the scan id column back into the data
dfCombattedCovPreserved['scan_id'] = scanIds
dfCombattedCovRemoved['scan_id'] = scanIds

# Save the combatted data both with and without covariate effects
dfCombattedCovPreserved.to_csv(phenoFn.replace('phenotypes_toCombat', 'combatted_covariates_preserved'), index=False)
dfCombattedCovRemoved.to_csv(phenoFn.replace('phenotypes_toCombat', 'combatted_covariates_removed'), index=False)

