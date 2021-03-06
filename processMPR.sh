#!/usr/bin/sh
#
# Usage:
# bash processMPR.sh file_list.txt

# Given a list of files
MPRLIST=$1
BASE=/cbica/projects/bgdimagecentral/22q11_dataorg

# Read the list of files
for fn in $(cat $MPRLIST) ; do
    tmp="${fn/rawdata/derivatives/freesurfer_6_0_0}"
    tmp2="${tmp//anat/}"
    FULLPATH=$(echo "$tmp2" | cut -f 1 -d '.')
    SUBJ=$(basename $FULLPATH)
    OUTDIR=$(dirname $FULLPATH)

    mkdir -p $OUTDIR

    # Preprocessing happens here
    echo $OUTDIR $SUBJ $fn

#    qsub $BASE/reconall-job.sh $OUTDIR $SUBJ $fn 
done
