%% Copy SDR files
% Carlos Rodriguez, Ph.D. Mind Research Network
% copies the output of the afni based ica_300_postproc script to compare to
% a previous version to perform correlations between voxel wise time
% courses.
% Run this script from the root testing directory

%% List the subjects directory
subjects = dir('sub-*');
subjects = subjects([subjects.isdir]); %modify subjects to contain only directories

%% For each subject, copy sdr files to afni_proc_files
for ii = 1:size(subjects,1)
    subj = subjects(ii).name
    copyfile([subj '/ses-baselineYear1Arm1/func/sdr*.nii'], 'afni_proc_files')
end

%% For each subject, delete the sdr files to run the icatb based regressed and detrended version of ica_300_postproc
for ii = 1:size(subjects,1)
    subj = subjects(ii).name
    delete([subj '/ses-baselineYear1Arm1/func/sdr*.nii'])
    delete([subj '/ses-baselineYear1Arm1/func/dr*.nii'])
    delete([subj '/ses-baselineYear1Arm1/func/fitts*.nii'])
    delete([subj '/ses-baselineYear1Arm1/func/r_sub*.nii'])
    delete([subj '/ses-baselineYear1Arm1/func/stats*'])
    delete([subj '/ses-baselineYear1Arm1/func/tcat*.nii'])
    delete([subj '/ses-baselineYear1Arm1/func/X_sub*'])
end
