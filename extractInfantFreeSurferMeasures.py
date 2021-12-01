import glob
import pandas as pd
import argparse
import os
from freesurfer_stats import CorticalParcellationStats


## 
# Get cortical measurements from both ?h.aparc.stats files
# @param fn A string representing the path to the original aseg file
# @param metric A string representing the metric name to get
# @return measure The numeric value of the specified measure or -1 if missing
def getMetricsFromAparcStats(fn, metric):
    # Get the two aparc.stats filenames from the base filename
    rhfn = fn.replace("brainvol.stats", "rh.aparc.stats")
    lhfn = fn.replace("brainvol.stats", "lh.aparc.stats")
    
    # Read the stats in for each hemisphere
    rhStats = CorticalParcellationStats.read(rhfn)
    lhStats = CorticalParcellationStats.read(lhfn)
    
    # Depending on the metric, need to look at different values
    if metric == "CSF":
        # Get whole brain values from the stats
        rhDf = rhStats.whole_brain_measurements
        lhDf = lhStats.whole_brain_measurements
        
        # Get left and right CSF
        rCsf = rhDf['brain_segmentation_volume_mm^3'] - rhDf['brain_segmentation_volume_without_ventricles_mm^3']
        lCsf = lhDf['brain_segmentation_volume_mm^3'] - lhDf['brain_segmentation_volume_without_ventricles_mm^3']
        
        # sum values together
        measure = rCsf + lCsf
        
    elif "Cortical" in metric:
        # Get all of the rows at the bottom of the aparg.stats file
        rhDf = rhStats.structural_measurements
        lhDf = lhStats.structural_measurements
        
        if "SurfaceArea" in metric:
            # Sum the values in the SurfArea column
            measure = sum(rhDf['surface_area_mm^2']) + sum(lhDf['surface_area_mm^2'])
        elif "ThickAvg":
            # Sum the values in the ThickAvg column
            measure = sum(rhDf['average_thickness_mm']) + sum(lhDf['average_thickness_mm'])
            
    else:
        measure = -1

    return measure


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
    # CLIP data has this structure
    if 'pat_id' in list(df):
        # Remove leading characters from subject id
        subjId = subj[4:]
        # Get the sex from the data frame
        sex = df[df['pat_id'] == subjId]['sex'].values[0]

    # Get the subject's sex
    elif 'subject_id' in list(df):
        print(subj)
        # Replace the subject id leading characters
        subjId = subj.replace('sub-22q', '22q_')
        # Get the sex from the data frame
        subjSubset = df[df['subject_id'] == subjId]
        if subjSubset.shape[0] > 0:
            sex = df[df['subject_id'] == subjId]['sex'].values[0]
        else:
            sex = 'U'

    return sex


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
    parser.add_argument('-f', '--fn', help='File containing demographic information about subject, specifically sex', required=True)

    args = parser.parse_args()

    path = args.dir
    demoFn = args.fn
    
    # Quick sanity check: does the input directory exist?
    if not os.path.exists(path):
        sys.exit("Error: the path doesn't exist:", path)

    # Quick sanity check: does the demographics file exist?
    if not os.path.exists(demoFn):
        sys.exit("Error: the demographics file doesn't exist:", demoFn)

    # After confirming the path does exist...

    # Initialize variables
    header = ["patient_id", "age_at_scan_days", "scan_id", "sex", "scanner_id", "BrainSeg", "CerebralWhiteMatter", "TotalGray", "EstimatedTotalIntraCranialVol", "SurfaceHoles", "SubCortGrayVol", "CSF", "SumCorticalSurfaceArea", "SumCorticalThickAvg"]
    fns = sorted(glob.glob(path+"/**/brainvol.stats", recursive=True))
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
        scannerId = str(scanId.split("FromScanner")[-1].split("_")[0])
        sex = getSexFromDemographics(demoDf, patId)
    
        # Read the aseg file
        with open(fn, 'r') as f:
            lines = f.readlines()
    
        # Create a new list
        row = [patId, patAge, scanId, sex, scannerId]
    
        # Pull out metrics we care about from the aseg file
        for metric in header[5:-6]:
            measure = getMeasureFromLine(lines, metric)
            row.append(measure)

        # Get two more measures not found in the brainvol.stats files
        row.append(getEstimatedTotalIntraCranialVol(fn))
        row.append(getEulerNumber(fn))

        # Add the subcortical gray volume
        measure = getMeasureFromLine(lines, header[-4])
        row.append(measure)

        # Add the three measures from aparc.stats
        for metric in header[-3:]:
            measure = getMetricsFromAparcStats(fn, metric)
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
