
#!/bin/bash

##################################################################################
# 03/22/2021
# Carlos Rodriguez, Ph.D., MRN

# The ica_300 postproc script utilizes two ways of extracting the motion Regressors
# from the .tsv files. The first way uses [trans_x..rot_z_derivative1_power2].
# The alternative uses, [trans_x..rot_z_power2], because the position of the last
# two columns can be switched depending on the subject.

# This script was created to double check the alternative method of selecting 
# regressor columns implemented in the ica_300 post fmriprep processing script.

# Usage: 1) navigate to the ica_300 log output folder 
#        2) bash check_regressors.sh
###################################################################################

# Make a list of the subject specific log files in the directory where the script 
# is called from
if [ ! -f tmp_subj_list.txt ]; then
    ls -f sub* > tmp_subj_list.txt
fi

# Create a file to track the subjects to check
touch tmp_check_regressors.txt

# Check each log file for the "Trying a different" string pattern displayed if the 
# alternative regressor extraction method was used, and out put to a temporary file
for i in `cat tmp_subj_list.txt`; do
    SUBJ=`echo ${i:0:19}`
    grep -q "Trying a different" $i
      if [ $? == 0 ]; then
        echo $SUBJ >> tmp_check_regressors.txt
      fi
done

# Count the number of line of subjects using the alt method
ALT=`wc -l < tmp_check_regressors.txt`
echo "Alternative regressor configuration was utilized in ${ALT} subjects."

#---------------------------------------------------------------------------------------
# For each subject in the check_regressors list, select the headers of the .tsv file,
# print out the column names of the extracted regressors, and count if any do not start
# with trans or rot

# Set Variables
SES="ses-baselineYear1Arm1"
#FILE="task-rest_run-1_desc-confounds_regressors.tsv"
DATAPATH="/export/research/analysis/human/jhouck/abcd/ica_300/ica_input"

touch tmp_tally.txt

# Go through each subject, and each run, and append the last
for SUBJ in `cat tmp_check_regressors.txt`; do    
    for FILE in `ls -f ${DATAPATH}/${SUBJ}/*regressors.tsv | cat`; do
        #echo ${FILE}
        TEST=`head -1 ${FILE} | tr $'\t' $'\n' | grep -A 20 trans_x | tail -1`
            if [[ $TEST == rot_z_power2 ]]; then
                echo ${TEST} >> tmp_tally.txt
            fi
    done
done

# Count the number of times rot_z_power2 was used for runs and divided by 4 for subjects
# will verify if each subject maintained a consistend confounds.tsv output structure and
# if the range from trans_x through rot_z_power2 was appropriate
COUNT_RUNS=`grep -c "rot_z_power2" tmp_tally.txt`
COUNT_SUB=`echo $COUNT_RUNS/4 | bc`

echo "${COUNT_RUNS} runs listed rot_z_power2 as the final column of the regressor extraction range."
echo "${COUNT_SUB} subjects listed rot_z_power2 as the final column of the regressor extraction range."

# Test if the number of subjects with alt log output equals the number of subjects with 
# rot_z_power2 as the last regressor in the range of 24 columns beginning with trans_x
if [[ $ALT == $COUNT_SUB ]]; then
    echo "Success! Alternative regressor extraction method verified."
fi

rm tmp*.txt

#-------------------------------------------------------------------------------------------
#In the fmriprep output, rot_z_derivative1_power2 can be swapped with rot_z_power2 as the last 
#column of the range of columns that begins with trans_x. The motion regressors are not in 
#fixed positions because acompcorr and steady state columns can vary. Despit this, all motition
#regressors are grouped together in a set of 24 adjacent columns, beginning with trans_x. 
#However, as noted, the last two columns can vary which complicates selecting the appropriate
#columns of motion regressors.
