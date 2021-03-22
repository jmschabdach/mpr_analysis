#!/usr/bin/sh
#
# Usage:
# bash processMPR.sh file_list.txt

# Given a list of files
MPRLIST=$1
BASE=/cbica/projects/bgdimagecentral/mpr_analysis

module load freesurfer/6.0.0

# Read the list of files
for fn in $(cat $MPRLIST) ; do

    # The next two lines assume input from /rawdata/sub/ses/anat
    tmp="${fn/rawdata/derivatives/mpr_preproc}"
    tmp2="${tmp/anat/}"

    FULLPATH="${tmp2/.nii.gz/}"

    echo "FULLPATH: $FULLPATH"
    echo ""

    mkdir -p $FULLPATH 

    # Preprocess the image
    qsub $BASE/preproc-job.sh $FULLPATH $fn $FULLPATH/preprocessed_output.nii.gz


done
