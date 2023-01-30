import os
import glob
import pandas as pd
import argparse

def parseSynthSegVolumes(ssFn):
    # load the subject's /synthseg_output/volumes.csv
    ssDf = pd.read_csv(ssFn)
     
#    ssDf['TCV'] = ssDf['total intracranial']
    ssDf['Cortex'] = ssDf['left cerebral cortex'] + ssDf['right cerebral cortex'] 
    ssDf['WMV'] = ssDf['left cerebral white matter'] + ssDf['right cerebral white matter']
    ssDf['TCV'] = ssDf['Cortex'] + ssDf['WMV']
    ssDf['sGMV'] = ssDf['left thalamus'] + ssDf['left caudate'] + ssDf['left putamen'] + ssDf['left pallidum'] + ssDf['left hippocampus'] + ssDf['left amygdala'] + ssDf['left accumbens area'] + ssDf['right thalamus'] + ssDf['right caudate'] + ssDf['right putamen'] + ssDf['right pallidum'] + ssDf['right hippocampus'] + ssDf['right amygdala'] + ssDf['right accumbens area']
#    ssDf['brainStem'] = ssDf['brain-stem']
    ssDf['Ventricles'] = ssDf['left lateral ventricle'] + ssDf['right lateral ventricle'] + ssDf['left inferior lateral ventricle'] + ssDf['right inferior lateral ventricle'] + ssDf['3rd ventricle'] + ssDf['4th ventricle'] 
    ssDf['CerebellumVolume'] = ssDf['left cerebellum cortex'] + ssDf['left cerebellum white matter'] + ssDf['right cerebellum cortex'] + ssDf['right cerebellum white matter']

    return ssDf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', help='Path to the directory containing the outputs of the SynthSeg 2.0 jobs')
    parser.add_argument('-f', '--demo-fn', help='Name of file containing demographics')

    args = parser.parse_args()
    print(args)

    path = args.directory
    outFn = os.path.join(path, "synthseg_2.0_phenotypes.csv")
    demoFn = args.demo_fn

    demoDf = pd.read_csv(demoFn)

    fsVersion = "SynthSeg_2.0"
    useFirstRunOnly = True
    phenoDfList = []

    for subj in sorted(os.listdir(path)):
        if "sub-" in subj:
            subjDir = os.path.join(path, subj)

            for sess in sorted(os.listdir(subjDir)):
                sessDir = os.path.join(subjDir, sess)
                # Pull out the age
                if "age" in sess:
                    ageDays = int(sess.split("age")[-1])
#                else:
#                    ageDays = "NA"
                # Get the list of volume output files
                volFns = sorted(glob.glob(os.path.join(sessDir, "**/volumes.csv"), recursive=True))
                print(subj, sess, ageDays, len(volFns))
                
#                if useFirstRunOnly:
#                    if type(volFns) is str:
#                        volFns = [volFns]
#                    else:
#                        volFns = [volFns[0]]
                for volFn in volFns:
                    print(volFn)
                    tmpDf = parseSynthSegVolumes(volFn)
                    tmpDf['subj_id'] = subj
                    tmpDf['sess_id'] = sess
                    tmpDf['scan_id'] = os.path.dirname(volFn)
                    tmpDf['fs_version'] = fsVersion
                    tmpDf['age_in_days'] = ageDays
                    print(tmpDf)
#
#                    # For CLIP
#                    s_id = subj[4:]
#                    if demoDf[demoDf['pat_id'] == s_id].shape[0] == 0:
#                        print(subj, "not in demographics file")
#                        reason = "MISSING"
#                        continue
#                    else:
#                        # Get the sex from the data frame
#                        tmpDf['sex'] = demoDf[demoDf['pat_id'] == s_id]['sex'].values[0]
#                        # Get the reason for the scan
#                        print(demoDf[demoDf['pat_id'] == s_id]['scan_reason_primary'].values[0])
#                        #print(demoDf[demoDf['pat_id'] == s_id]['scan_reason_categories'].values)
#                        #print(demoDf[demoDf['pat_id'] == s_id]['confirm_neurofibromatosis'].values)
#
#                        tmpDf['scan_reason_primary'] = demoDf[demoDf['pat_id'] == s_id]['scan_reason_primary'].values[0]
#                        print(tmpDf['scan_reason_primary'])
#                        tmpDf['scan_reason_categories'] = sorted(demoDf[demoDf['pat_id'] == s_id]['scan_reason_categories'].values)[-1]
#                        tmpDf['confirm_neurofibromatosis'] = sorted(demoDf[demoDf['pat_id'] == s_id]['confirm_neurofibromatosis'].values)[-1]
#                        # Get the scan grade
#                        tmpDf['rawdata_image_grade'] = demoDf[demoDf['pat_id'] == s_id]['rawdata_image_grade'].values[0]
#                    
#                        # Add demographic info
#                        if '1p5' in s_id:
#                            tmpDf['MagneticFieldStrength'] = 1.5
#                        elif '3p0' in s_id:
#                            tmpDf['MagneticFieldStrength'] = 3.0
#
#                        tmpDf['scanner_id'] = [str(os.path.dirname(volFn).split("FromScanner")[-1].split("_")[0]) for volFn in volFns]
#
#                        # Rearrange columns
#                        cols = list(tmpDf)
#                        cols = cols[-12:] + cols[:-12]
#                        tmpDf = tmpDf[cols]
#                        print(tmpDf)
#                        # Add the dataframe for a single scan to the list of dataframes
#                        phenoDfList.append(tmpDf)
#    
#    groupDf = pd.concat(phenoDfList, ignore_index=True)
#    groupDf.to_csv(outFn, index=False)


if __name__ == "__main__":
    main()
