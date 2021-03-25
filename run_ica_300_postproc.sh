#!/bin/bash

##############################################################################
#This script was created to run the ica_300_postproc matlab function in a loop 
#to process subjects resting state fmri data. The script relies on matlab and
#AFNI programs. The post processing involves determining the number of rest 
#files, selecting the files with the least amount of motion, gathering motion
# regressors from the confounds.tsv files, performing motion regression of 
# voxelwise time series, despiking, and smoothing.

# Input argument is the root directory of the fmriprepped abcd data
# Usage: bash run_ica_300_postproc.sh /export/research/analysis/human/jhouck/abcd/ica_300/fmriprep_out/fmriprep
 
##############################################################################
RDR=$1

# Change directory into the root directory containing all fmriprepped subjects
cd ${RDR}

# list the sub* directories, pipe through sed to remove the trailing slash 
# then redirect to directories.txt

ls -d sub*/ | sed 's/.$//' > /export/research/analysis/human/jhouck/abcd/ica_300/ica_input/directories.txt

#begin for loop to process each of the sample subjects
for SUBJ in `cat /export/research/analysis/human/jhouck/abcd/ica_300/ica_input/directories.txt` ; do
	echo $SUBJ
	SUBJ_DIR=$RDR/$SUBJ
	matlab -nodisplay -nojvm -r "afni_ica_300_postproc('$SUBJ_DIR'); exit"
done

#remove the directories.txt file

#rm directories.txt
