#!/bin/sh
#
# The name of the job
#$ -N synthseg
#
# Join stdout and stderrL:
#$ -j y
#
# Set the amount of memory being requested.
#$ -l h_vmem=16G

FN=$1

source ~/.software/miniconda3/etc/profile.d/conda.sh

# Prep the environment
module unload freesurfer/5.3

export FREESURFER_HOME=/cbica/projects/bgdimagecentral/.software/freesurfer-dev-2022-04-11
source $FREESURFER_HOME/SetUpFreeSurfer.sh

conda activate synthseg


# Prep the filename

OUT=${FN/rawdata/derivatives/mpr_fs_synthsegplus}
OUT=${OUT/.nii.gz//synthseg_volume.nii.gz}
METRICS=${OUT/synthseg_volume.nii.gz/metrics}


echo "mri_synthseg --i $FN --o $OUT --vol $METRICS --robust"

mri_synthseg --i $FN --o $OUT --vol $METRICS --robust
