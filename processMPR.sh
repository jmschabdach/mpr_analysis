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

    tmp="${fn/rawdata/derivatives/freesurfer_6_0_0_preproc}"
    tmp2="${tmp//anat/}"
    echo $tmp2
    FULLPATH="${tmp2/.nii.gz/}"
    SUBJ=$(basename $FULLPATH)
    OUTDIR=$(dirname $FULLPATH)

    mkdir -p $FULLPATH 

    # Preprocess the image
    bash preprocWashUACPCAlignment.sh --workingdir=$OUTDIR/preprocessing --in=$fn --out=$OUTDIR/preprocessing/preprocessed_output.nii.gz --omat="premat.mat"

    # Run recon-all
#    qsub $BASE/reconall-job.sh $OUTDIR $SUBJ $fn 

done
