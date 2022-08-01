from neuroHarmonize import harmonizationLearn
import pandas as pd
import numpy as np

# Load the data
phenoFn = '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_phenotypes_toCombat.csv'
covFn = '/Users/youngjm/Data/clip/fs6_stats/fs6_structural_stats_covariates_toCombat.csv'
# Old data, sanity check on code running
#phenoFn = "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc0-2_phenotypes_toCombat.csv"
#covFn = "/Users/youngjm/Data/clip/images/derivatives/mpr_fs_reconall_6.0.0_tables/mpr_fs_reconall_6.0.0_clip_qc0-2_covariates_toCombat.csv"

phenoDf = pd.read_csv(phenoFn)
covDf = pd.read_csv(covFn)
print(phenoDf.shape)
# debugging: take the first 10 rows
phenoDf = phenoDf.drop(index=[379])
covDf = covDf.drop(index=[379])
#phenoDf = phenoDf.iloc[:380, ]
#covDf = covDf.iloc[:380, ]

scanIds = phenoDf['scan_id']
data = np.array(phenoDf.drop(columns=['scan_id']))
#data = np.array(phenoDf[['TotalGrayVol', 'eTIV']])

# Specify categorical columns
catCols = ['sex', 'reason']

# Set up the covariates
covDf = pd.get_dummies(covDf, columns=catCols, drop_first=True)

# Run ComBat
print(type(data))
print(data.shape)
print(type(covDf))
print(covDf.shape)

with np.errstate(invalid='ignore', divide='ignore'):
    model, data_combatted, mystery = harmonizationLearn(data, covDf, return_s_data=True)

print(type(model))
print(data_combatted.shape)
print(data_combatted[0])
print(mystery.shape)
print(mystery[0])

