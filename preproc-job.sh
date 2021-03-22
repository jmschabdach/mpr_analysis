#!/bin/sh
#
# The name of the job
#$ -N preprocess-mpr-job
#
# Join stdout and stderrL:
#$ -j y
#
# Specify environment variables
# -v SUBJECTS_DIR="/cbica/projects/bgdimagecentral/images/derivatives/freesurfer_7_1_0/"
# -v FREESURFER_HOME="/cbica/software/external/freesurfer/centos7/7.1.0/"

#module purge
module unload freesurfer/5.3.0
module load freesurfer/6.0.0

# bash /cbica/software/external/freesurfer/centos7/7.1.0/SetUpFreeSurfer.sh

source $FREESURFER_HOME/SetUpFreeSurfer.sh

BASE=/cbica/projects/bgdimagecentral/mpr_analysis

WORKINGDIR=$1
IN=$2
OUT=$3

bash $BASE/preprocWashUACPCAlignment.sh --workingdir=$WORKINGDIR --in=$IN --out=$OUT --omat="premat.mat"

