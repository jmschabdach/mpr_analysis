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


## Load the number of holes from the surf/lh.orig.euler.txt and surf/rh.orig.euler.txt
# @param fn A string specifying the brainvol.stats file location
# @return Sum of the number of holes in the lh and rh surface files
def getEulerNumber(fn):
    eulerPath = os.path.dirname(os.path.dirname(fn))
    eulerPath = os.path.join(eulerPath, "surf")
    lhEulerFn = os.path.join(eulerPath, "lh.orig.euler.txt")
    rhEulerFn = os.path.join(eulerPath, "rh.orig.euler.txt")

    # Open the files and load the lines
    with open(lhEulerFn, 'r') as f:
        lhLines = f.readlines()

    with open(rhEulerFn, 'r') as f:
        rhLines = f.readlines()

    # Get the first line in the file containing "holes"
    lhMatchingLine = [i for i in lhLines if "holes" in i][0]
    rhMatchingLine = [i for i in rhLines if "holes" in i][0]

    # Extract the number of holes from each line
    lhEulerNum = int(lhMatchingLine.strip().split("-->")[-1].split("holes")[0].strip())
    rhEulerNum = int(rhMatchingLine.strip().split("-->")[-1].split("holes")[0].strip())

    return lhEulerNum + rhEulerNum


## Load the stats/lh.aparc.stats file and get the eTIV value from it
# @param fn A string specifying the brainvol.stats file
# @returns eTIV The measure parsed from the file
def getEstimatedTotalIntraCranialVol(fn):
    # variable set up: the file is in the same directory as the brainvol.stats file
    statsPath = os.path.dirname(fn)
    lhStatsFn = os.path.join(statsPath, "lh.aparc.stats")

    # Open the file and read the contents
    with open(lhStatsFn, 'r') as f:
        lhLines = f.readlines()

    # Parse the measure from the file contents
    eTIV = getMeasureFromLine(lhLines, "EstimatedTotalIntraCranialVol")

    return eTIV


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', help='Directory containing the outputs of a FreeSurfer reconall pipeline', required=True)

    args = parser.parse_args()

    path = args.dir
    
    # Quick sanity check: does the input directory exist?
    if not os.path.exists(path):
        sys.exit("Error: the path doesn't exist:", path)

    # After confirming the path does exist...

    # Initialize variables
    header = ["patient_id", "age_at_scan_days", "scan_id", "BrainSeg", "CerebralWhiteMatter", "TotalGray", "EstimatedTotalIntraCranialVol", "SurfaceHoles"]
    fns = sorted(glob.glob(path+"/**/brainvol.stats", recursive=True))
    newRows = []

    # Quick sanity check: does the input directory contain aseg.stats files somewhere?
    if not len(fns) > 0:
        sys.exit("Error: the directory does not contain aseg.stats files")

    # for each file
    for fn in fns:
        # get info that's important for consideration in analyses
        patId = fn.split("/")[-3].split("_")[0]
        patAge = str(int(fn.split("/")[-3].split("_")[1].split("age")[-1]))
        scanId = fn.split("/")[-3]
    
        # Read the aseg file
        with open(fn, 'r') as f:
            lines = f.readlines()
    
        # Create a new list
        row = [patId, patAge, scanId]
    
        # Pull out metrics we care about from the aseg file
        for metric in header[3:-2]:
            measure = getMeasureFromLine(lines, metric)
            row.append(measure)

        # Get two more measures not found in the brainvol.stats files
        row.append(getEstimatedTotalIntraCranialVol(fn))
        row.append(getEulerNumber(fn))
    
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
