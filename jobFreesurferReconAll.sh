#!/bin/sh
#
# The name of the job
#$ -N reconall
#
# Join stdout and stderrL:
#$ -j y
#
# Set the amount of memory being requested.
#$ -l h_vmem=8G

#source $FREESURFER_HOME/SetUpFreeSurfer.sh

OUTDIR=$1     # string specifying where the output of recon-all should live when it's done
INPUT=$2      # file name
SUBJ=$3       # identifier string for the subject (sub-XXXX)
FSLVERSION=$4 # string

BASE=/cbica/projects/bgdimagecentral/Projects/mpr_analysis_jenna
SUBJECTS_DIR=$OUTDIR

source /cbica/projects/bgdimagecentral/.software/miniconda3/etc/profile.d/conda.sh

echo "Inputs:"
echo "OUTDIR: $OUTDIR"
echo "INPUT: $INPUT"
echo "SUBJ: $SUBJ"
echo "FSLVERSION: $FSLVERSION"

if [[ $FSLVERSION == "6.0.0" ]] ; then
    module unload freesurfer/5.3.0
    module load freesurfer/6.0.0
elif [[ $FSLVERSION == "5.3.0" ]] ; then
    module unload freesurfer/6.0.0
    module load freesurfer/5.3.0
elif [[ $FSLVERSION == "7.1.1" ]] ; then
    module load freesurfer/7.1.0
else
    echo "Something wrong with FSLVERSION $FSLVERSION"
fi

echo "Scratch directory: "

# Run the preprocessing
SCRATCH=${OUTDIR/reconall/preprocWashUACPCAlignment}
echo "SCRATCH: $SCRATCH"

INTERIM=$SCRATCH/preprocessed_mpr.nii.gz
# do I need to copy the input image into the working directory?
mkdir -p $SCRATCH
cp $INPUT $SCRATCH
# Check orientation
fslreorient2std $INPUT $SCRATCH/orig_mni152_aligned.nii.gz
# Deface 
mri_deface $SCRATCH/orig_mni152_aligned.nii.gz /cbica/projects/bgdimagecentral/templates/talairach_mixed_with_skull.gca /cbica/projects/bgdimagecentral/templates/face.gca $SCRATCH/defaced.nii.gz
# Use the HCP preprocessing step (yes that's actually where it's from)
bash $BASE/preprocWashUACPCAlignment.sh --workingdir=$SCRATCH --in=$SCRATCH/defaced.nii.gz --out=$INTERIM --omat="premat.mat"

# Run recon-all 
recon-all -sd $OUTDIR -subject $SUBJ -i $INTERIM -all -target /cbica/projects/bgdimagecentral/.fsaverage
#recon-all -sd $OUTDIR -subject $SUBJ -i $INPUT -all -target /cbica/projects/bgdimagecentral/.fsaverage

 
# Job ended, clean up
rm $OUTDIR/fsaverage
mv $OUTDIR/$SUBJ/* $OUTDIR
rm -r $OUTDIR/$SUBJ
