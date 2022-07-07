# Usage: submit a series of FS and IFS reconall jobs
# bash [filename] [BIDS rawdata directory] [integer freesurfer version]

DATADIR=$1
FS=$2         # 6 or 7

# Specify the FS version to load on CUBIC
if [ $FS -eq 6 ] ; then
    FSVERSION="6.0.0"
elif [ $FS -eq 7 ] ; then
    FSVERSION="7.1.0"
fi

BASE=/cbica/projects/bgdimagecentral/Projects/mpr_analysis_jenna

# for subject in directory
for subj in $DATADIR/sub-* ; do
    # for session in directory
    SUBJID=$(basename $subj)
    for session in $subj/ses-* ; do

        echo $subj

        # Get the age of the subject at the time of the scan
        sesIdStr=${session#*ses-}
        ageStr=${sesIdStr#*age}
        ageNum=$((10#$ageStr))


        # Check that an MPR exists in the session 
        for fn in $session/anat/*.nii.gz ; do 
            if [[ ${fn,,} == *"mpr"* ]] ; then      # Cool an MPR exists
           
                # if the age is less than or equal to 3 years
                if [ $ageNum -lt 1096 ] ; then
                    # Set up the output directory 
                    OUTDIR="${fn/rawdata/derivatives/mpr_ifs_reconall_defaced_$FSVERSION}"
                    OUTDIR="${OUTDIR/anat/}"
                    OUTDIR=${OUTDIR%%.nii.gz}
                    echo $OUTDIR

                    # If the output directory doesn't exist, then make it
                    if [ ! -d $OUTDIR ] ; then
                        mkdir -p $OUTDIR
                    fi

                    # Submit the IFS job
#                    qsub $BASE/jobInfantFreesurferReconAll.sh $fn $OUTDIR $ageNum $FSVERSION

                # else age > 3 years
                else 
                    # Set up the output directory
                    OUTDIR="${fn/rawdata/derivatives/mpr_fs_reconall_defaced_$FSVERSION}"
                    OUTDIR="${OUTDIR/anat/}"
                    OUTDIR=${OUTDIR%%.nii.gz}
                    echo $OUTDIR
       
                    # If the ouput directory doesn't exist, then make it
                    if [ ! -d $OUTDIR ] ; then
                        mkdir -p $OUTDIR
                    fi

                    # Submit the FS job
                    qsub $BASE/jobFreesurferReconAll.sh $OUTDIR $fn $SUBJID $FSVERSION
                fi
            fi
        done # end for fn in session/anat/*.nii.gz
    done # end for session in subject
done # end for subject
