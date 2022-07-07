#!/bin/sh
#
# The name of the job
#$ -N reconall
#
# Join stdout and stderrL:
#$ -j y
#
# Set the amount of memory being requested.
#$ -l h_vmem=4G
#
# Specify environment variables
# -v SUBJECTS_DIR="/cbica/projects/bgdimagecentral/images/derivatives/freesurfer_7_1_0/"
# -v FREESURFER_HOME="/cbica/software/external/freesurfer/centos7/7.1.0/"

#module purge
# bash /cbica/software/external/freesurfer/centos7/7.1.0/SetUpFreeSurfer.sh

source $FREESURFER_HOME/SetUpFreeSurfer.sh

OUTDIR=$1
SUBJ=$2
INPUT=$3
FSLVERSION=$4

if [[ $FSLVERSION == "6.0.0" ]] ; then
    module unload freesurfer/5.3.0
    module load freesurfer/6.0.0
elif [[ $FSLVERSION == "5.3.0" ]] ; then
    module load freesurfer/5.3.0
fi

echo "Scratch directory: "
echo $SBIA_TMPDIR

# Move files to tmp dir
SUBJECTS_DIR=$SBIA_TMPDIR

# mkdir $SUBJECTS_DIR/$SUBJ
cp $INPUT $SUBJECTS_DIR
cp -r /cbica/projects/bgdimagecentral/fsaverage/ $SUBJECTS_DIR

# Set a trap to copy any temp files 
run_on_exit(){
    cp -r $SUBJECTS_DIR/$SUBJ/* $OUTDIR
}
trap run_on_exit EXIT

echo "Contents of directory pre-run: "
ls $SBIA_TMPDIR
 
recon-all -subject $SUBJ -i $INPUT -all -target $SUBJECTS_DIR/fsaverage

cp -r $SUBJECTS_DIR/$SUBJ/* $OUTDIR

#recon-all -subject sub-22q0002 -i /cbica/projects/bgdimagecentral/images/rawdata/sub-22q0002/ses-11957/anat/sub-22q0002_ses-11957_acq-MPR09Iso_run-001_T1w.nii.gz -all
# /cbica/software/external/freesurfer/centos7/6.0.0/bin/recon-all -subject sub-22q0002 -i /cbica/projects/bgdimagecentral/images/rawdata/sub-22q0002/ses-11957/anat/sub-22q0002_ses-11957_acq-MPR09Iso_run-001_T1w.nii.gz -all
