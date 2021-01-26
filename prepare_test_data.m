%% Prepare test data
% Carlos Rodriguez, Ph.D. Mind Research Network
% Randomly select 10 subjects to copy into a testing directory for the
% ica_300_postproc function.

% Change directory to the fmri prepped subject directories
cd('/export/research/analysis/human/jhouck/abcd/ica_300/fmriprep_out/fmriprep')

% List the subject folders
subjects = dir('sub-*');
subjects = subjects([subjects.isdir]); %modify subjects to contain only directories

%randomly select 10 subjects
sel_sub = randi([1, size(subjects,1)], 1,10);

%copy the directories and the .html files to ica_300_postproc testing
%directory
for ii = 1:size(sel_sub,2)
    
    if isfile(['/export/research/analysis/human/jhouck/abcd/ica_300/fmriprep_out/fmriprep/' subjects(sel_sub(ii)).name '.html'])
    copyfile(['/export/research/analysis/human/jhouck/abcd/ica_300/fmriprep_out/fmriprep/' subjects(sel_sub(ii)).name '.html'], '/export/research/analysis/human/jhouck/cobre06_65007/carlos_work/ica_300/ica_300_postproc')
    copyfile(['/export/research/analysis/human/jhouck/abcd/ica_300/fmriprep_out/fmriprep/' subjects(sel_sub(ii)).name], ['/export/research/analysis/human/jhouck/cobre06_65007/carlos_work/ica_300/ica_300_postproc/' subjects(sel_sub(ii)).name])

    end
end

% Move to the testing directory to view resuls
cd('/export/research/analysis/human/jhouck/cobre06_65007/carlos_work/ica_300/ica_300_postproc')


