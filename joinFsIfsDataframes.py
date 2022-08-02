import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--ifs-fn', help='Path to InfantFreeSurfer output file', required=True)
    parser.add_argument('-f', '--fs-fn', help='Path to FreeSurfer output file', required=True)
    parser.add_argument('-o', '--outfn', help='Path to desired output file')

    args = parser.parse_args()

    fsFn = args.fs_fn
    ifsFn = args.ifs_fn
    outFn = args.outfn

    # Load the tables
    fsDf = pd.read_csv(fsFn)
    ifsDf = pd.read_csv(ifsFn)
    print(fsDf.shape)
    print(ifsDf.shape)

    # Drop duplicates
    fsDf = fsDf.drop_duplicates()
    ifsDf = ifsDf.drop_duplicates()
    print(fsDf.shape)
    print(ifsDf.shape)

    # Add column for processing info
    fsDf['fs_version'] = 'FS_6.0.0'
    ifsDf['fs_version'] = 'IFS_6.0.0'

    # Make the confirm_neurofibromatosis column consistent
    ifsDf['confirm_neurofibromatosis'] = ifsDf['confirm_neurofibromatosis'].fillna(0.0)
    fsDf['confirm_neurofibromatosis'] = fsDf['confirm_neurofibromatosis'].fillna(0.0)

    # Print the number of columns that are not shared
    print(len([i for i in list(fsDf) if i not in list(ifsDf)]))
    print([j for j in list(ifsDf) if j not in list(fsDf)])

    # Combine the dataframes
    newDf = pd.concat([fsDf, ifsDf], axis=0)
    print(newDf.shape)
    
    # Drop any columns with missing values (ie in only one df)
    newDf = newDf.dropna(axis=1)

    # Save the new dataframe
    newDf.to_csv(outFn)


if __name__ == "__main__":
    main()

