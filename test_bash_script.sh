#!/bin/bash

#This script was created to test the implementation of the
#ica_300_postproc matlab function in a loop to process test
#subjects hoping to catch any errors that may arise.

# list the sub* directories, pipe through sed to remove the trailing slash then redirect to directories.txt
ls -d sub*/ | sed 's/.$//' > directories.txt

# begin for loop to process each of the sample subjects
for subj in `cat ./directories.txt` ; do
	echo $subj
	matlab -nodisplay -nojvm -r "ica_300_postproc('/export/research/analysis/human/jhouck/cobre06_65007/carlos_work/ica_300/ica_300_postproc/$subj');exit"
done

#remove the directories.txt file
rm directories.txt
