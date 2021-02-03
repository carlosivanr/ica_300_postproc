#!/bin/bash

#This script was created to correlate the output of the
#ica_300_postproc matlab function that uses afni commands only
#with the icatb_ica_300_postproc matlab function that uses the
#regress and detrend functions from icatb and the old afni 3d
#despike.

# list the sub* directories, pipe through sed to remove the trailing slash then redirect to directories.txt
ls -d sub*/ | sed 's/.$//' > temp_directories.txt

# begin for loop to process each of the sample subjects
for subj in `cat ./temp_directories.txt` ; do
	echo $subj
	matlab -nodisplay -nojvm -r "icatb_ica_300_postproc('/export/research/analysis/human/jhouck/cobre06_65007/carlos_work/ica_300/ica_300_postproc/$subj');exit"
	matlab -nodisplay -nojvm -r "time_series_correlations('/export/research/analysis/human/jhouck/cobre06_65007/carlos_work/ica_300/ica_300_postproc/$subj');exit"
done

#remove the directories.txt file
rm temp_directories.txt
