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

def extractAparcPhenotypes(fn, side):
    columnHeaders = []
    values = []

    # load the file using the regular file open
    fn = fn.replace('aseg', side+'.aparc')
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
        columnHeaders.append(columnName)
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
        columnHeaders.append(feature+"_surfArea")
        values.append(surfArea)
        columnHeaders.append(feature+"_grayVol")
        values.append(grayVol)
        columnHeaders.append(feature+"_thickAvg")
        values.append(thickAvg)

    return columnHeaders, values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', help='Directory containing the outputs of a FreeSurfer reconall pipeline', required=True)

    args = parser.parse_args()

    path = args.dir
    
    # Quick sanity check: does the input directory exist?
    if not os.path.exists(path):
        sys.exit("Error: the path doesn't exist:", path)

    fns = sorted(glob.glob(path+"/**/aseg.stats", recursive=True))

    # Quick sanity check: does the input directory contain aseg.stats files somewhere?
    if not len(fns) > 0:
        sys.exit("Error: the directory does not contain aseg.stats files")

    # create a blank dataframe
    mainDf = pd.DataFrame()

    # for each file
    for fn in fns:
        # get info that's important for consideration in analyses
        patId = fn.split("/")[-3].split("_")[0]
        scanId = fn.split("/")[-3]

        # Get the aseg.stats phenotypes
        asegHeaders, asegValues = extractAsegPhenotypes(fn)

        # Get the aparc.stats phenotypes
        lhHeaders, lhValues = extractAparcPhenotypes(fn, 'lh')
        rhHeaders, rhValues = extractAparcPhenotypes(fn, 'rh')

        # Make the scan dataframe
        scanDf = pd.DataFrame(asegValues+lhValues+rhValues, keys=asegHeaders+lhHeaders+rhHeaders)

        # Add the scan dataframe to the main dataframe
        mainDf = pd.concat([mainDf, scanDf], ignore_index=True)

    
    # Concatenate all of the rows into a dataframe
    df = pd.DataFrame(newRows, columns=header)
     
    # Save the dataframe
    outFn = os.path.join(os.path.dirname(path), os.path.basename(path)+"_structural_stats.csv")

    df.to_csv(outFn, index=False)

    # Let the user know the data has been extracted and saved
    print("The data from", path, "has been saved to", outFn)


if __name__ == "__main__":
    main()
