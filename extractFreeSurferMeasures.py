import glob
import pandas as pd
import argparse
import os
import sys
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
def getMetricsFromAparcStats(fn, metric):
    print("Getting metrics from aparc stats - this is a flag that this function is ever called")
    # Get the two aparc.stats filenames from the base filename
    rhfn = fn.replace("aseg.stats", "rh.aparc.stats")
    lhfn = fn.replace("aseg.stats", "lh.aparc.stats") 

    # Read the stats in for each hemisphere
    rhStats = CorticalParcellationStats.read(rhfn)
    lhStats = CorticalParcellationStats.read(lhfn)

    # Depending on the metric, need to look at different values
    if metric == "CSF":
        # Get whole brain values from the stats
        rhDf = rhStats.whole_brain_measurements
        lhDf = lhStats.whole_brain_measurements

        # Get left and right CSF
        rCsf = rhDf['brain_segmentation_volume_mm^3'].values[0] - rhDf['brain_segmentation_volume_without_ventricles_mm^3'].values[0]
        lCsf = lhDf['brain_segmentation_volume_mm^3'].values[0] - lhDf['brain_segmentation_volume_without_ventricles_mm^3'].values[0]

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

            if metric.startswith("Avg"):
                denom = len(rhDf['average_thickness_mm']) + len(lhDf['average_thickness_mm'])
                measure = measure/denom

    else:
        measure = -1

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', help='Directory containing the outputs of a FreeSurfer reconall pipeline', required=True)
    parser.add_argument('-i', '--infant', help='Flag to indicate that the directory contains Infant FreeSurfer outputs', action='store_true')

    args = parser.parse_args()

    path = args.dir
    isIfs = args.infant

    baseFn = "aseg.stats"
    if isIfs:
        baseFn = "brainvol.stats"
    
    # Initialize variables
    fns = sorted(glob.glob(path+"/**/"+baseFn, recursive=True))
    print(fns)
    newRows = []

    # Quick sanity check: does the input directory contain aseg.stats files somewhere?
    if not len(fns) > 0:
        sys.exit("Error: the directory does not contain "+baseFn +" files")

    # Create a blanket dataframe
    mainDf = pd.DataFrame()

    demoHeader = ['subjId', 'scanId']

    # for each file
    for fn in fns:
        # get info that's important for consideration in analyses
        scanId = os.path.dirname(fn.replace('/stats/','/' )).split("/")[-1]
        if fn.endswith("/"):
            subjId = fn.replace(path, '').split("/")[0]
        else:
            subjId = fn.replace(path, '').split("/")[1]
        print(path)
        print(fn)
        print(scanId)
        print(subjId)
 
        #THIS IS A CRAPPY TEMPORARY FIX BUT IT WORKS
        if 'sub-' in scanId:
            subjId = scanId.split('_')[0].split('-')[-1]

        demoValues = [subjId, scanId]
        
        # Read the aseg file
        with open(fn, 'r') as f:
            lines = f.readlines()
     
        # Get the aseg.stats phenotypes
        asegHeaders, asegValues = extractAsegPhenotypes(fn)

        # Get the aparc.stats phenotypes
        lhHeaders, lhValues = extractAparcPhenotypes(fn, 'lh', baseFn)
        rhHeaders, rhValues = extractAparcPhenotypes(fn, 'rh', baseFn)
        
        headers = demoHeader+asegHeaders+lhHeaders+rhHeaders
        values = demoValues+asegValues+lhValues+rhValues

        # Get a few more parameters that are missing in IFS
        if isIfs:
            etiv = getEstimatedTotalIntraCranialVolIfs(fn)
            surfaceHoles = getEulerNumberIfs(fn)
            rhMeanThickness = getCorticalThicknessIfs(fn, 'rh', baseFn)
            lhMeanThickness = getCorticalThicknessIfs(fn, 'lh', baseFn)
            headers = demoHeader+asegHeaders+lhHeaders+rhHeaders+['eTIV', 'SurfaceHoles', 'rh_MeanThickness', 'lh_MeanThickness']
            values = demoValues+asegValues+lhValues+rhValues+[etiv, surfaceHoles, rhMeanThickness, lhMeanThickness]

        # Make the scan dataframe
        dataDict = {}
        for key, value in zip(headers, values):
            dataDict[key] = [value]

        scanDf = pd.DataFrame(data=dataDict)

        # Add the scan dataframe to the main dataframe
        mainDf = pd.concat([mainDf, scanDf], ignore_index=True)


    # Add missing columns that have been key in other analyses
    if "rh_BrainSegVol" in list(mainDf) and "rh_BrainSegVolNotVent" in list(mainDf):
        mainDf['VentricleVolume'] = (mainDf['rh_BrainSegVol'].astype(float) - mainDf['rh_BrainSegVolNotVent'].astype(float)) + (mainDf['lh_BrainSegVol'].astype(float) - mainDf['lh_BrainSegVolNotVent'].astype(float))

    if "rh_WhiteSurfArea" in list(mainDf):
        mainDf['CorticalSurfaceArea'] = mainDf['rh_WhiteSurfArea'].astype(float) + mainDf['lh_WhiteSurfArea'].astype(float) #???

    if "rh_MeanThickness" in list(mainDf):
        mainDf['MeanCorticalThickness'] = (mainDf['rh_MeanThickness'].astype(float) + mainDf['lh_MeanThickness'].astype(float))/2
                
    # Save the dataframe
    outFn = os.path.join(os.path.dirname(path), os.path.basename(path)+"_structural_stats.csv")

    mainDf.to_csv(outFn, index=False)

    # Let the user know the data has been extracted and saved
    print("The data from", path, "has been saved to", outFn)


if __name__ == "__main__":
    main()
