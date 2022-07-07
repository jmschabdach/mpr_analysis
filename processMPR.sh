#!/usr/bin/sh
#
# Usage:
# bash processMPR.sh file_list.txt

# Given a list of files
MPRLIST=$1
FSLVERSION=$2
BASE=/cbica/projects/bgdimagecentral/Projects/mpr_analysis

if [[ $FSLVERSION == "6.0.0" ]] ; then
    module unload freesurfer/5.3.0
    module load freesurfer/6.0.0
elif [[ $FSLVERSION == "5.3.0" ]] ; then
    module load freesurfer/5.3.0
fi

# Read the list of files
for fn in $(cat $MPRLIST) ; do

    tmp=$fn

    if [[ "$tmp" == *"anat"* ]] ; then
        tmp="${tmp//anat/}"
    fi

    # This section executes if the input has been preprocessed/lives in the preprocessed folder
    if [[ "$tmp" == *"mpr_preproc"* ]] ; then
        if [[ $FSLVERSION == "5.3.0" ]] ; then
            if [[ $fn == *"1p5"* ]] ; then
                tmp="${tmp/rawdata/derivatives/mpr_post_reconall_1p5T_FSL5p3}"
            elif [[ $fn == *"3p0"* ]] ; then
                tmp="${tmp/rawdata/derivatives/mpr_post_reconall_3p0T_FSL5p3}"
            fi
        elif [[ $FSLVERSION == "6.0.0" ]] ; then
            if [[ $fn == *"1p5"* ]] ; then
                tmp="${tmp/rawdata/derivatives/mpr_post_reconall_1p5T_FSL6p0}"
            elif [[ $fn == *"3p0"* ]] ; then
                tmp="${tmp/rawdata/derivatives/mpr_post_reconall_3p0T_FSL6p0}"
            fi
        fi

        FULLPATH=$(dirname $tmp)
        OUTDIR=$(dirname $FULLPATH)
        SUBJ=$(basename $(dirname $OUTDIR))
 

        echo "file:     $tmp"
        echo "SUBJ:     $SUBJ"
        echo "OUTDIR:   $OUTDIR"
        echo "FULLPATH: $FULLPATH"
    
        mkdir -p $FULLPATH 
    
        qsub $BASE/reconall-job.sh $FULLPATH $SUBJ $fn $FSLVERSION

    else 
        echo "Nope, data not preprocessed: $fn"
    fi

done
