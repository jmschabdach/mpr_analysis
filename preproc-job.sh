#!/bin/sh
#
# The name of the job
#$ -N preprocess-mpr-job
#
# Join stdout and stderrL:
#$ -j y
#
# Specify environment variables
## -v SUBJECTS_DIR="/cbica/projects/bgdimagecentral/images/derivatives/freesurfer_7_1_0/"
## -v FREESURFER_HOME="/cbica/software/external/freesurfer/centos7/7.1.0/"

#module purge
#module unload freesurfer/5.3.0
#module load freesurfer/6.0.0

#source $FREESURFER_HOME/SetUpFreeSurfer.sh

BASE=/cbica/projects/bgdimagecentral/Projects/mpr_analysis

WORKINGDIR=$1
IN=$2
OUT=$3
FSLVERSION=$4


if [[ $FSLVERSION == "5.3.0" ]] ; then
    module unload freesurfer/6.0.0
    module load freesurfer/5.3.0
elif [[ $FSLVERSION == "6.0.0" ]] ; then
    module unload freesurfer/5.3.0
    module load freesurfer/6.0.0
fi

bash $BASE/preprocWashUACPCAlignment.sh --workingdir=$WORKINGDIR --in=$IN --out=$OUT --omat="premat.mat"

