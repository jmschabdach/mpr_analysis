#!/bin/sh
#
# The name of the job
#$ -N ifs-reconall
#
# Join stdout and stderrL:
#$ -j y
#
# Set the amount of memory being requested
#$ -l h_vmem=32G

INPUT=$1
OUTDIR=$2
AGE=$3
FSL_VERSION=$4

# Set up the environment
module unload freesurfer/5.3.0


# Set up environment variable and directory structure IFS is expecting
SUBJECTS_DIR=$(dirname $OUTDIR)
SUBJ=$(basename $OUTDIR)
echo $SUBJECTS_DIR
ls $SUBJECTS_DIR
ls $SUBJECTS_DIR/$SUBJ

# Copy the input image to the output directory
cp $INPUT $OUTDIR/mprage.nii.gz

# Convert age in days to age in months
AGE_MONTHS=$( awk -v var1=$AGE -v daysPerYear="365.25" -v monthsPerYear="12" 'BEGIN { print ( int(( var1 / daysPerYear ) * monthsPerYear ) ) }')
echo "Age in months: $AGE_MONTHS"

# Do stuff to deal with the different FS versions
if [[ $FSL_VERSION == *"7.1.0"* ]] ; then
    # Set up the Freesurfer environment variables
    export FREESURFER_HOME=/cbica/projects/bgdimagecentral/.software/freesurfer-7.1.1/
    INFANT_FREESURFER=/cbica/projects/bgdimagecentral/.software/freesurfer-infant-7

elif [[ $FSL_VERSION == *"6.0.0"* ]] ; then
    #module load freesurfer/6.0.0
    #export FREESURFER_HOME=/cbica/projects/bgdimagecentral/.software/freesurfer-6.0.0/
    export FREESURFER_HOME=/cbica/projects/bgdimagecentral/.software/freesurfer-infant-6
    INFANT_FREESURFER=/cbica/projects/bgdimagecentral/.software/freesurfer-infant-6
fi

bash $INFANT_FREESURFER/SetUpFreeSurfer.sh
source $INFANT_FREESURFER/FreeSurferEnv.sh
PATH="$PATH:$FREESURFER_HOME:/cbica/projects/bgdimagecentral/.software/miniconda/bin/perl"
export PATH
export SUBJECTS_DIR

# Run recon-all
$INFANT_FREESURFER/bin/infant_recon_all --s $SUBJ --age $AGE_MONTHS --stats --ccseg

# Remove the mprage.nii.gz file from the output directory
# rm $OUTDIR/mprage.nii.gz

