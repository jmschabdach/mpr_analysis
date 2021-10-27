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

    tmp=$fn

    # The next two lines assume input from /rawdata/sub/ses/anat
    if [[ "$tmp" == *"rawdata"* ]] ; then
        tmp="${tmp/rawdata/derivatives/freesurfer_6_0_0_preproc}"
    elif [[ "$tmp" == *"mpr_preproc"* ]] ; then
        tmp="${tmp/mpr_preproc/freesurfer_6_0_0_preproc}"
    fi 

    if [[ "$tmp" == *"anat"* ]] ; then
        tmp="${tmp//anat/}"
    fi

#    FULLPATH="${tmp/.nii.gz/}"
#    SUBJ=$(basename $FULLPATH)
#    OUTDIR=$(dirname $FULLPATH)

    FULLPATH=$(dirname $tmp)
    OUTDIR=$(dirname $FULLPATH)
    SUBJ=$(basename $(dirname $OUTDIR))

    echo "SUBJ:     $SUBJ"
    echo "OUTDIR:   $OUTDIR"
    echo "FULLPATH: $FULLPATH"

    mkdir -p $FULLPATH 

    qsub $BASE/reconall-job.sh $FULLPATH $SUBJ $fn


done
