from neuroHarmonize import harmonizationLearn
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--phenotypes')
    parser.add_argument('-c', '--covariates')
    parser.add_argument('-o', '--output-fn')

    args = parser.parse_args()

    # Load the data
    # phenoFn = '/Users/youngjm/Data/chop-nf1/fs6_stats/fs6_structural_stats_phenotypes_toCombat.csv'
    # covFn = '/Users/youngjm/Data/chop-nf1/fs6_stats/fs6_structural_stats_covariates_toCombat.csv'
    phenoFn = args.phenotypes
    covFn = args.covariates
    outFn = args.output_fn

    phenoDf = pd.read_csv(phenoFn) 
    covDf = pd.read_csv(covFn)
    print(phenoDf.shape)
    print(covDf.shape)
    
    # REMOVE THIS FOR OTHER DATA, ISSUE WITH A SINGLE DATASET (CLIP)
    #x=379
    #print(phenoDf.iloc[x,])
    #print(covDf.iloc[x,])
    #phenoDf = phenoDf.drop(index=[x])
    #covDf = covDf.drop(index=[x])
    #phenoDf = phenoDf.iloc[x:, ]
    #covDf = covDf.iloc[x:, ]

    #phenoDf = phenoDf[['scan_id', 'MeanCorticalThickness']]
    
    print(phenoDf.shape)
    print(covDf.shape)
    
    # Pull out the scan ids from the dataframe
    scanIds = phenoDf['scan_id']
    data = np.array(phenoDf.drop(columns=['scan_id']))
    
    # Specify categorical columns
    catCols = ['sex', 'reason']
    
    # Set up the covariates
    covDf = pd.get_dummies(covDf, columns=catCols, drop_first=True)
#    np.seterr(invalid='ignore')
    
    # Run ComBat
    # -  Using the return_s_data flag returns a third arg: combatted data with covariate effects not preserved
    model, data_combatted = harmonizationLearn(data, covDf) #, smooth_terms=['log_age_in_days'])
    
    # Convert the data (numpy.array) into dataframes (pd.DataFrame)
    dfCombattedCovPreserved = pd.DataFrame(columns=list(phenoDf.drop(columns=['scan_id'])), data=data_combatted)
    
    # Because we're potentially ignoring a division error, add a line to drop anything in the combatted data with nan
#    print(dfCombattedCovPreserved.shape)
#    dfCombattedCovPreserved = dfCombattedCovPreserved.dropna(axis=0)
#    print(dfCombattedCovPreserved.shape)
    
    # Add the scan id column back into the data
    dfCombattedCovPreserved['scan_id'] = scanIds
    
    # Save the combatted data both with and without covariate effects
    dfCombattedCovPreserved.to_csv(outFn, index=False)

if __name__ == "__main__":
    main()
