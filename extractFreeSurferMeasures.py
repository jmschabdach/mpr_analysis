import glob
import pandas as pd
import argparse
import os


## Extract a single measure from the top of an aseg.stats file
# @param lines A list of lines read in from a *.stats file
# @param measureName A string specifying the name of the measure to extract
# @returns measure The extracted measure as a string
def getMeasureFromLine(lines, measureName):
    # Find the line containing the measure name
    matchingLine = [l for l in lines if measureName in l]

    # Parse the measure out of the matching line
    measure = matchingLine[0].split(",")[-2]

    # Return the extracted measure
    return measure


## Get the sex of a subject from the demographics file
# @param df A dataframe containing the demographic information
# @param subj A string specifying the subject id
# @returns sex The sex of the subject as a string
def getSexFromDemographics(df, subj):
    # Remove leading characters from subject id
    subjId = subj[4:]
    # Get the subject's sex
    sex = df[df['pat_id'] == subjId]['sex'].values[0]

    return sex


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', help='Directory containing the outputs of a FreeSurfer reconall pipeline', required=True)
    parser.add_argument('-f', '--fn', help='File containing demographic information about subject, specifically sex', required=True)

    args = parser.parse_args()

    path = args.dir
    demoFn = args.fn
    
    # Quick sanity check: does the input directory exist?
    if not os.path.exists(path):
        sys.exit("Error: the path doesn't exist:", path)

    # Quick sanity check: does the demographic file exist?
    if not os.path.exists(demoFn):
        sys.exit("Error: the demographics file doesn't exist:", demoFn)

    # After confirming the path does exist...

    # Initialize variables
    header = ["patient_id", "age_at_scan_days", "sex", "scan_id", "BrainSeg", "CerebralWhiteMatter", "TotalGray", "EstimatedTotalIntraCranialVol", "SurfaceHoles"]
    fns = sorted(glob.glob(path+"/**/aseg.stats", recursive=True))
    newRows = []

    # Load the demographics file
    demoDf = pd.read_csv(demoFn)

    # Quick sanity check: does the input directory contain aseg.stats files somewhere?
    if not len(fns) > 0:
        sys.exit("Error: the directory does not contain aseg.stats files")

    # for each file
    for fn in fns:
        # get info that's important for consideration in analyses
        patId = fn.split("/")[-3].split("_")[0]
        patAge = str(int(fn.split("/")[-3].split("_")[1].split("age")[-1]))
        scanId = fn.split("/")[-3]
        sex = getSexFromDemographics(demoDf, patId)
    
        # Read the aseg file
        with open(fn, 'r') as f:
            lines = f.readlines()
    
        # Create a new list
        row = [patId, patAge, scanId, sex]
    
        # Pull out metrics we care about from the aseg file
        for metric in header[4:]:
            measure = getMeasureFromLine(lines, metric)
            row.append(measure)
    
        # Add new row to list of rows
        newRows.append(row)
    
    # Concatenate all of the rows into a dataframe
    df = pd.DataFrame(newRows, columns=header)
     
    # Save the dataframe
    outFn = os.path.join(os.path.dirname(path), os.path.basename(path)+"_structural_stats.csv")

    df.to_csv(outFn, index=False)

    # Let the user know the data has been extracted and saved
    print("The data from", path, "has been saved to", outFn)


if __name__ == "__main__":
    main()
