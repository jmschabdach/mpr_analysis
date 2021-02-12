#!/usr/bin/sh
#
# Usage:
# bash processMPR.sh file_list.txt

# Given a list of files
MPRLIST=$1
BASE=/cbica/projects/bgdimagecentral/mpr_analysis

# Read the list of files
for fn in $(cat $MPRLIST) ; do
    tmp="${fn/rawdata/derivatives/freesurfer_6_0_0}"
    tmp2="${tmp//anat/}"
    FULLPATH=$(echo "$tmp2" | cut -f 1 -d '.')
    SUBJ=$(basename $FULLPATH)
    OUTDIR=$(dirname $FULLPATH)

    mkdir -p $FULLPATH 


    echo $FULLPATH
    echo $SUBJ $fn  

    qsub $BASE/reconall-job.sh $FULLPATH $SUBJ $fn 
done
