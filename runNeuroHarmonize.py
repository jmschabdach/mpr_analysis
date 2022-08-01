from neuroHarmonize import harmonizationLearn
import pandas as pd
import numpy as np

# Load the data
phenoFn = '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_phenotypes_toCombat.csv'
covFn = '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_covariates_toCombat.csv'

phenoDf = pd.read_csv(phenoFn)
covDf = pd.read_csv(covFn)

scanIds = phenoDf['scan_id']
data = np.array(phenoDf.drop(columns=['scan_id']))

# Specify categorical columns
catCols = ['sex', 'reason']

# Set up the covariates
covDf = pd.get_dummies(covDf, columns=catCols, drop_first=True)

# Run ComBat
print(type(data))
print(data.shape)
print(type(covDf))
print(covDf.shape)
model, data_combatted, mystery = harmonizationLearn(data, covDf, return_s_data=True)

print(type(model))
print(data_combatted.shape)
print(mystery.shape)
