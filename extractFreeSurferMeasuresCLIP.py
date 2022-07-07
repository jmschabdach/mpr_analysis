import glob
import pandas as pd
import argparse
import os
from freesurfer_stats import CorticalParcellationStats


def extractAsegPhenotypes(fn):
    columnHeaders = []
    values = []

    # load the file using the regular file open
    with open(fn, 'r') as f:
        lines = f.readlines()

    # for any line starting with "# Measure"
    measureLines = [row for row in lines if "# Measure" in row]
    for row in measureLines:
        # split by ", "
        rowEls = row.split(", ")
        columnName = rowEls[1]
        value = rowEls[-2]

        # Add the stuff to the lists
        columnHeaders.append(columnName)
        values.append(value)
        
    # get the local volumetric measures
    measureLines = [row for row in lines if "# " not in row]
    for row in measureLines:
        # split on whitespace
        rowEls = row.split()
        columnName = rowEls[4]
        value = rowEls[3]
        
        # Add the stuff to the lists
        columnHeaders.append(columnName)
        values.append(value)

    return columnHeaders, values


def extractAparcPhenotypes(fn, side, statsBase):
    columnHeaders = []
    values = []

    # load the file using the regular file open
    fn = fn.replace(statsBase, side+'.aparc.stats')
    with open(fn, 'r') as f:
        lines = f.readlines()

    # for any line starting with "# Measure Cortex"
    measureLines = [row for row in lines if "# Measure" in row]
    for row in measureLines:
        # split by ", "
        rowEls = row.split(", ")
        columnName = rowEls[1]
        value = rowEls[-2]

        # Add the stuff to the lists
        columnHeaders.append(side+"_"+columnName)
        values.append(value)
        
    # get the local volumetric measures
    measureLines = [row for row in lines if "# " not in row]
    for row in measureLines:
        # split on whitespace
        rowEls = row.split()
        feature = rowEls[0]
        surfArea = rowEls[2]
        grayVol = rowEls[3]
        thickAvg = rowEls[4]
        
        # Add the stuff to the lists
        columnHeaders.append(side+"_"+feature+"_surfArea")
        values.append(surfArea)
        columnHeaders.append(side+"_"+feature+"_grayVol")
        values.append(grayVol)
        columnHeaders.append(side+"_"+feature+"_thickAvg")
        values.append(thickAvg)

    return columnHeaders, values

##
# Get cortical measurements from both ?h.aparc.stats files
# @param fn A string representing the path to the original aseg file
# @param metric A string representing the metric name to get
# @return measure The numeric value of the specified measure or -1 if missing
def getCorticalThicknessIfs(fn, side, statsBase):
    # Get the two aparc.stats filenames from the base filename
    newFn = fn.replace(statsBase, side+".aparc.stats")

    # Read the stats in for each hemisphere
    stats = CorticalParcellationStats.read(newFn)

    # Get all of the rows at the bottom of the aparg.stats file
    df = stats.structural_measurements

    #      if "SurfaceArea" in metric:
    #          # Sum the values in the SurfArea column
    #          measure = sum(rhDf['surface_area_mm^2']) + sum(lhDf['surface_area_mm^2'])
    # Sum the values in the ThickAvg column
    measure = sum(df['average_thickness_mm'])
    denom = len(df['average_thickness_mm'])
    measure = measure/denom

    return measure

## Extract a single measure from the top of an aseg.stats file
# @param lines A list of lines read in from a *.stats file
# @param measureName A string specifying the name of the measure to extract
# @returns measure The extracted measure as a string
def getMeasureFromLine(lines, measureName, matchIdx=0):
    # Find the line containing the measure name
    matchingLine = [l for l in lines if measureName in l]

    # Parse the measure out of the matching line
    measure = matchingLine[matchIdx].split(",")[-2]

    # Return the extracted measure
    return measure


## Load the number of holes from the surf/lh.orig.euler.txt and surf/rh.orig.euler.txt
# @param fn A string specifying the brainvol.stats file location
# @return Sum of the number of holes in the lh and rh surface files
def getEulerNumberIfs(fn):
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
def getEstimatedTotalIntraCranialVolIfs(fn):
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
    parser.add_argument('-m', '--main', help='File containing demographics information, specifically primary scan reason', required=True)
    parser.add_argument('-i', '--infant', help='Flag indicating the directory is for Infant FreeSurfer outputs', action='store_true')

    args = parser.parse_args()

    path = args.dir
    mainFn = args.main
    isIfs = args.infant
    
    # Quick sanity check: does the input directory exist?
    if not os.path.exists(path):
        sys.exit("Error: the path doesn't exist:", path)

    # Quick sanity check: does the demographic file exist?
    if not os.path.exists(mainFn):
        sys.exit("Error: the demographics file doesn't exist:", mainFn)

    # After confirming the path does exist...

    # Initialize variables
    demoHeaders = ["patient_id", "age_at_scan_days", "scan_id", "sex", "MagneticFieldStrength", "scanner_id", "scan_reason_primary", "scan_reason_categories", "confirm_neurofibromatosis", "rawdata_image_grade"]# "BrainSeg", "CerebralWhiteMatter", "TotalGray", "EstimatedTotalIntraCranialVol", "SurfaceHoles", "SubCortGrayVol", "CSF", "Cortex", "SumCorticalSurfaceArea", "SumCorticalThickAvg", "AvgCorticalThickAvg"]

    statsBase = "aseg.stats"
    if isIfs:
        statsBase = "brainvol.stats"

    fns = sorted(glob.glob(path+"/**/"+statsBase, recursive=True))
    newRows = []

    # Load the demographics file
    masterDf = pd.read_csv(mainFn)

    # Quick sanity check: does the input directory contain aseg.stats files somewhere?
    if not len(fns) > 0:
        sys.exit("Error: the directory does not contain aseg.stats files")

    # Create a blanket dataframe
    mainDf = pd.DataFrame()

    # for each file
    for fn in fns:
        # get info that's important for consideration in analyses
        patId = fn.split("/")[-3].split("_")[0]
        patAge = str(int(fn.split("/")[-3].split("_")[1].split("age")[-1]))
        scanId = fn.split("/")[-3]
        if '1p5' in scanId:
            fieldStrength=1.5
        elif '3p0' in scanId:
            fieldStrength=3.0
        scannerId = str(scanId.split("FromScanner")[-1].split("_")[0])
        
        # Remove leading characters from subject id
        subjId = patId[4:]

        if masterDf[masterDf['pat_id'] == subjId].shape[0] == 0:
            print(subjId, "not in demographics file")
            reason = "MISSING"
            continue

        else:
            # Get the sex from the data frame
            sex = masterDf[masterDf['pat_id'] == subjId]['sex'].values[0]
            # Get the reason for the scan
            reason = masterDf[masterDf['pat_id'] == subjId]['scan_reason_primary'].values[0]
            reasons = masterDf[masterDf['pat_id'] == subjId]['scan_reason_categories'].values[0]
            nf = masterDf[masterDf['pat_id'] == subjId]['confirm_neurofibromatosis'].values[0]
            # Get the scan grade
            grade = masterDf[masterDf['pat_id'] == subjId]['rawdata_image_grade'].values[0]
    
        # Read the aseg file
        with open(fn, 'r') as f:
            lines = f.readlines()
    
        # Create a new list
        demoValues = [patId, patAge, scanId, sex, fieldStrength, scannerId, reason, reasons, nf, grade]
    
        # Get the aseg.stats phenotypes
        asegHeaders, asegValues = extractAsegPhenotypes(fn)

        # Get the aparc.stats phenotypes
        lhHeaders, lhValues = extractAparcPhenotypes(fn, 'lh', statsBase)
        rhHeaders, rhValues = extractAparcPhenotypes(fn, 'rh', statsBase)

        headers = demoHeaders+asegHeaders+lhHeaders+rhHeaders
        values = demoValues+asegValues+lhValues+rhValues

        # Get a few more parameters that are missing in IFS
        if isIfs:
            etiv = getEstimatedTotalIntraCranialVolIfs(fn)
            surfaceHoles = getEulerNumberIfs(fn)
            rhMeanThickness = getCorticalThicknessIfs(fn, 'rh', statsBase)
            lhMeanThickness = getCorticalThicknessIfs(fn, 'lh', statsBase)
            headers = demoHeaders+asegHeaders+lhHeaders+rhHeaders+['eTIV', 'SurfaceHoles', 'rh_MeanThickness', 'lh_MeanThickness']
            values = demoValues+asegValues+lhValues+rhValues+[etiv, surfaceHoles, rhMeanThickness, lhMeanThickness]

        # Make the scan dataframe
        dataDict = {}
        for key, value in zip(headers, values):
            dataDict[key] = [value]

        scanDf = pd.DataFrame(data=dataDict)

        # Add the scan dataframe to the main dataframe
        mainDf = pd.concat([mainDf, scanDf], ignore_index=True)


    # Add missing columns that have been key in other analyses
    if 'rh_BrainSegVol' in list(mainDf):
        mainDf['VentricleVolume'] = (mainDf['rh_BrainSegVol'].astype(float) - mainDf['rh_BrainSegVolNotVent'].astype(float)) + (mainDf['lh_BrainSegVol'].astype(float) - mainDf['lh_BrainSegVolNotVent'].astype(float))
   
    if 'rh_WhiteSurfArea' in list(mainDf):
        mainDf['CorticalSurfaceArea'] = mainDf['rh_WhiteSurfArea'].astype(float) + mainDf['lh_WhiteSurfArea'].astype(float) #???


    if 'rh_MeanThickness' in list(mainDf):
        mainDf['MeanCorticalThickness'] = (mainDf['rh_MeanThickness'].astype(float) + mainDf['lh_MeanThickness'].astype(float))/2.0


    # Save the dataframe
    outFn = os.path.join(os.path.dirname(path), os.path.basename(path)+"_structural_stats.csv")

    mainDf.to_csv(outFn, index=False)

    # Let the user know the data has been extracted and saved
    print("The data from", path, "has been saved to", outFn)


if __name__ == "__main__":
    main()
