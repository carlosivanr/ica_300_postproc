# ica_300_postproc
Matlab script that processes fmriprep output of ABCD data

The collection of scripts in this repository were created to process resting state (rs) fMRI data from the ABCD study. The ABCD study is the largest longitudinal study of adolescent brain developement, enrolling close to 12,000 participants to be tracked for a period of 10 years. 

The scripts are Matlab and shell based to apply post processing to the rs-fMRI data had been pre processed with fMRI Prep. The post processing steps include:
- verifying files
- despiking
- smoothing

Other scripts were created to correlate the output of the ica_300_postproc matlab function that uses afni commands only with the icatb_ica_300_postproc matlab function that uses the regress and detrend functions from icatb and the old afni 3d despike function to verify pipelines.
