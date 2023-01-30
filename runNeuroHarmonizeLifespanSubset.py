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
    phenoFn = args.phenotypes
    covFn = args.covariates
    outFn = args.output_fn

    phenoDf = pd.read_csv(phenoFn) 
    covDf = pd.read_csv(covFn)
    print(phenoDf.shape)
    print(covDf.shape)

    phenoDf = phenoDf.astype({'scan_id':'string'})
    covDf = covDf.astype({'SITE': 'string'})
    
    # REMOVE THIS FOR OTHER DATA, ISSUE WITH A SINGLE DATASET (CLIP)
#    toDrop = [i for i in range(50, 100)] + [i for i in range(200, 400)] 
#    toDrop += [i for i in range(450, 500)]
#    x=600
#    y=700
#    #print(phenoDf.iloc[x,])
#    #print(covDf.iloc[x,])
#    #phenoDf = phenoDf.drop(index=[x])
#    #covDf = covDf.drop(index=[x])
#    phenoDf = phenoDf.iloc[x:y, ]
#    covDf = covDf.iloc[x:y, ]

#    dupCovRows = covDf.duplicated()
#    idxToDrop = [i for i, x in enumerate(dupCovRows) if x]
#    covDf = covDf.drop(idxToDrop, axis=0)
#    phenoDf = phenoDf.drop(idxToDrop, axis=0)
#    scanIds = scanIds.drop(idxToDrop, axis=0)
#    print("Dropping rows with repeat covariates")
#
#    dupPhenoRows = phenoDf.drop(columns=['scan_id']).duplicated()
#    idxToDrop = [i for i, x in enumerate(dupPhenoRows) if x]
#    covDf = covDf.drop(idxToDrop, axis=0)
#    phenoDf = phenoDf.drop(idxToDrop, axis=0)
#    scanIds = scanIds.drop(idxToDrop, axis=0)
#    print("Dropping rows with repeat phenotypes")


    #phenoDf = phenoDf[['scan_id', 'MeanCorticalThickness']]
    
    print(phenoDf.shape)
    print(covDf.shape)
    
    # Pull out the scan ids from the dataframe
    scanIds = phenoDf['scan_id']
    data = np.array(phenoDf.drop(columns=['scan_id']))

    # Specify categorical columns
    catCols = ['sex', 'fs_version']
    
    # Set up the covariates
    covDf = pd.get_dummies(covDf, columns=catCols, drop_first=True)
##    np.seterr(invalid='ignore')
    
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
    dfCombattedCovPreserved['scan_id'] = scanIds.values
    print(dfCombattedCovPreserved.head())

    # Save the combatted data both with and without covariate effects
    dfCombattedCovPreserved.to_csv(outFn, index=False)

if __name__ == "__main__":
    main()
