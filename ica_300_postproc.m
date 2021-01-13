function [subj_dir_path] = ica_300_postproc(subj_dir_path)

%% ICA 300 post fmriprep processing function
% Carlos Rodriguez, Ph.D. 
% This function was created for the ica_300 project to unzip the best of 4 
% runs of resting state fmri from 300 participants of the ABCD study.
% Voxelwise timecourses are detrended using gift icatb_detren
% Motion parameters are then regressed out of the voxelwise timecourses
% using gift icatb_regress. 
%
% Dependencies include: AFNI version 20.2.02
% Set AFNI_AUTOGZIP = NO in .afnirc. 
% Matlab parallel pool for despiking and smoothing multiple subjects at once.
% Regression is handled in serial.
%
%% Usage from matlab command line:
% argument in = absolute path to the subject directory
% ex.
% ica_300_postproc('/export/reasearch/analysis/human/jhouck/abcd/ica_300/sub-NDAR*')
% requires write permissions for error_logging in the directory from which
% it is called.
%
%% Usage from bash terminal:
% matlab -nodisplay -nojvm -r "ica_300_postproc('absolute/path/to/subject-specific-fmriprep-output-directory');exit;"
%
%% Revision History
% 12/28/20 - Updated 3dTcat and parallel processing for smoothing
% 12/29/20 - Removed 3dTcat for concatenating, added code to read .tsv 
%   confound regressors. Added detrend and regress functions, file merge of
%   3d volumes, 3dDespike
% 12/30/20 - Added display messages, cleaned up file merge code section. Modified
% code to run from subjects 2 through 8 in testing.
%   check timepoints after regression
%   check slices in 3dDespike
% 01/04/21 - Modified reading of .tsv files to use readtable instead of
% tdfread which speeds up the process, but still requires conversion of
% 'n/a's to double. Set AFNI 3dDespike prefix to include .nii to bypass the
% .brik to .nii conversion step. Added -nomask and -ssave options. 
% Smoothing set to run within each subject folder.
% 01/05/21 - Modified script to work as a matlab function that takes the subject 
% specific directory as an argument. Changed preproc_bold to
% preproc_brain_bold.nii which are masked files
% 01/06/21 - Added parfor loop functionality to 3dDespike, added
% functionality to use fslmerge if FSLOUTPUTTYPE is set to NIFTI.
% 01/07/21 - Added option to smooth with AFNI
% 01/08/21 - Added code to delete pre-existing smoothed files to prevent
% AFNI errors.
% 01/11/21 - Modified file and fmriprep error checks. Added code to set
% environmental variables for AFNI, OMP_NUM_THREADS, and FSLOUTPUTTYPE.
% Script will now look at time points to exclude runs from further
% processing.
% 01/12/21 - Replaced icatb and SPM functions with AFNI routines. Added 
% code to unzip mask files to use in AFNI routines.

% Possible modifications to make:
% - add session input argument, in case there will be session 2 data to
% process
% - add error log output argument for directing where to save files
% - parfor loop for unzipping files

%% Code section for testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subjects = dir('sub-NDAR*'); %list all of the files and folders that begin with sub-* prefix
% subjects = subjects([subjects.isdir]); %modify list to contain only directories
% subj_dir_path  = [subjects(7).folder filesep subjects(7).name]; %set subj_dir_path to subject

%% Set environmental variables
tic
setenv('PATH', [getenv('PATH') ':/export/research/analysis/human/jhouck/shared/tools/abin_20202']); %Set AFNI Version by appending to the path
setenv('OMP_NUM_THREADS', '16'); %Set # of openmp threads

%% Verify directories, files, and errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ses = 'ses-baselineYear1Arm1';
[parent_dir_path, subj, ~] = fileparts(subj_dir_path); %retrieves parent directory and subject ID
disp(['Processing ' subj])
disp('Checking func directory ...')
% check if the func directory exists
if isfolder([subj_dir_path filesep ses filesep 'func'])
    disp('func directory found. Checking files...')
    
    % check the number of rs files
    rs_runs = dir([subj_dir_path filesep ses '/func/*preproc_bold_brain.nii.gz']); %list all the zipped rest runs
    masks = dir([subj_dir_path filesep ses '/func/*brain_mask.nii.gz']);
    if size(rs_runs, 1) >= 4
        disp('At least 4 resting state files found. Checking fmriprep html report for errors...')
        
        % check the html report for errors
        if isfile([parent_dir_path filesep subj '.html']) %if file exists, proceed to the next lines
            %T.exists_html(i) = 1; %fill in the table with the .html exists
            [status, cmdout] = system(['grep "No errors to report!" ' parent_dir_path filesep subj '.html']); % search the html report
            
            %status of grep for "No errors found!", 0=no error/successful, 1=not found, 2=permission denied for grep            
            if contains(cmdout, 'No errors') %status == 0
                disp('No errors found. Proceeding to post fmriprep processing.')
            elseif status == 1 %grep for "No errors" returns no results, check for errors
                [~, cmdout] = system(['grep "Errors" ' parent_dir_path filesep subj '.html']); % search the html report
                if contains(cmdout, 'Errors') %
                    disp('Errors in fmriprep html report found. Exiting ica_300_postproc.')
                    log_msg = [subj ', Errors in fmriprep html report found'];
                    fid = fopen('post_proc_log_file.txt', 'a'); %open the log file to append message
                    fprintf(fid, '%s: %s\n', datestr(now, 0), log_msg);
                    fclose(fid);
                    return
                end
            elseif status == 2
                disp('func directory exists and at least 4 rs-fMRI files found, but fmriprep html report is unverifiable. Possible permission error.')
                log_msg = [subj ', fmriprep html report unverifiable, possible permission error'];
                fid = fopen('post_proc_log_file.txt', 'a'); %open the log file to append message
                fprintf(fid, '%s: %s\n', datestr(now, 0), log_msg);
                fclose(fid);      
            end
        else
            disp('No fmriprep html report found, errors unverifiable. Exiting ica_300_postproc.')
                log_msg = [subj ', fmriprep html report not found'];
                fid = fopen('post_proc_log_file.txt', 'a'); %open the log file to append message
                fprintf(fid, '%s: %s\n', datestr(now, 0), log_msg);
                fclose(fid);
           return
        end       
    else
        %disp([ num2str(size(rs_runs,1)) ' resting state files found. Exiting ica_300_postproc.'])
        fprintf(2, [ num2str(size(rs_runs,1)) ' resting state files found. Exiting ica_300_postproc.\n']) %prints error message in red
        log_msg = [subj ', '  num2str(size(rs_runs,1)) ' resting state files found.'];
        fid = fopen('post_proc_log_file.txt', 'a');
        fprintf(fid, '%s: %s\n', datestr(now, 0), log_msg); %open the log file to append message
        fclose(fid);
        return
        
    end
else
    disp('No func directory found. Exiting ica_300_postproc.')
    log_msg = [subj ', No func directory found.'  ];
    fid = fopen('post_proc_log_file.txt', 'a');
    fprintf(fid, '%s: %s\n', datestr(now, 0), log_msg); %open the log file to append message
    fclose(fid);
    return
end
clear status

%% Conditional statement to proceed if func directory exists, at least 4 rs-fMRI runs exist, and no error reported by fmriprep
if contains(cmdout, 'No errors')
    % Change directory to ica_300 fmriprep output directory
    cd([subj_dir_path filesep ses filesep 'func']);%

    % Testing mode, clears files
    disp('Testing mode: Clearing files before processing.')
    delete sdr_*.nii; delete dr_*.nii; delete r_*.nii; delete tcat_*.nii; delete stats_*.nii; delete fitts_*.nii; delete X_sub*.xmat.1D;
    
    % QC number of runs, timepoints, and number of motion outliers
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(rs_runs,1) >= 4
        disp('Checking which files to process.')
        tsv_files = dir('*_regressors.tsv'); %list .tsv confound regressor files
        num_outliers = zeros(1, size(tsv_files, 1)); %preallocate number of motion outliers per run
        num_timepoints = zeros(1, size(tsv_files, 1)); %preallocate number of time points per run
        for ii = 1:size(tsv_files,1)
            [~, mot] = system(['grep -o "motion_outlier" ' tsv_files(ii).name ' | wc -l']);
            num_outliers(ii) = str2double(mot);
            [~, ntp] = system(['wc -l <' tsv_files(ii).name]);
            num_timepoints(ii) = str2double(ntp);
        end
        run_info = [1:size(tsv_files,1); num_timepoints; num_outliers]'; % Create a matrix with the run_info
        T = array2table(run_info);  % Write run info to a .csv file as part of logging information regarding file processing
        T.Properties.VariableNames = {'run', 'timepoints', 'num_outliers'};
        writetable(T, [subj '_run_info.csv']); %writes out a table with run diagnostics
        disp(T); clear T
        run_info(run_info(:, 2) ~= 383, :) = [] ; %deletes any rows with timepoints not equal to 383.
    end
    
    % If remaining runs with 383 timepoints are less than 4, exit out of
    % the function and write to error log
    if size(run_info, 1) < 4
        disp('Not enough runs to process. Number of timepoints may be less than expected. Exiting ica_300_postproc.')
        log_msg = [subj ', Not of enough runs to process. Number of timepoints may be less than expected.'];
        fid = fopen('post_proc_log_file.txt', 'a');
        fprintf(fid, '%s: %s\n', datestr(now, 0), log_msg); %open the log file to append message
        fclose(fid);
        return
    elseif size(run_info, 1) >= 4 %rank the remaining files according to the number of motion outliers
        run_info = sortrows(run_info, 3); %sort run_info according to the number of motion outliers
        runs_to_use = run_info(1:4, 1); %select the the first 4 runs which contain the least number of motion outliers
        runs_to_use = sort(runs_to_use); %sort runs_to_use to reflect run acquisition
        disp(['Using runs ' num2str(runs_to_use(1)) ', ' num2str(runs_to_use(2)) ', ' num2str(runs_to_use(3)) ', and ' num2str(runs_to_use(4)) '.'])
    end
    clear run_info num_time_points num_outliers

    % Unzip resting state files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:4 %this loop uses runs_to_use as an index to unzip files
        filezp = rs_runs(runs_to_use(jj)).name;
        fileun = filezp(1:end-3); %
        run = rs_runs(runs_to_use(jj)).name(strfind(rs_runs(runs_to_use(jj)).name, 'run'):strfind(rs_runs(runs_to_use(jj)).name, 'run')+4); %identifies which run
        if isfile(fileun)
            disp([run ' rest file unzipped.'])
        else
           disp([run ' rest file is zipped. Unzipping...']) 
           system(['gunzip -k ' filezp]);
        end
    end
    
    % Unzip mask files
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:4 %this loop uses runs_to_use as an index to unzip files
        filezp = masks(runs_to_use(jj)).name;
        fileun = filezp(1:end-3); %
        run = ['run-' num2str(runs_to_use(jj))];
        if isfile(fileun)
            disp([run ' mask file unzipped.'])
        else
           disp([run ' mask file is zipped. Unzipping...']) 
           system(['gunzip -k ' filezp]);
        end
    end     
    
    % Regression of motion (performed in serial)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Performing regression of voxelwise timecourses.')
    
    % Initialize empty cells for file checking
    reg_files = cell(4,1);
    dr_files = cell(4,1);
    sdr_files = cell(4,1);
    mask_files = cell(4,1);
    
    tic; 
    for kk = 1:4
        run = ['run-' num2str(runs_to_use(kk))]; %identifies which run is processing
        file_name = [subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %file name to detrend & regress
        mask_file = [subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-brain_mask.nii'];%
        tsv_file_name = [subj '_' ses '_task-rest_' run '_desc-confounds_regressors.tsv']; %concatenates a .tsv file name
        
        reg_files{kk,1} = ['r_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for despiking regressed files
        dr_files{kk,1} = ['dr_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for smoothing despiked files
        sdr_files{kk,1} = ['sdr_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for smoothing despiked files
        mask_files{kk,1} = mask_file; %initialize a cell for despiking

        % check if files exists
        if isfile(reg_files{kk})
            disp([run ' has been regressed. Skipping...'])
        else
            disp([run ' has not been regressed. Extracting regressors and performing regression of motion parameters.'])
            
            % Motion Regression with 3dDeconvolve
            subj_run = [subj '_' run];
            
            % Prepare a 1d file with the motion regressors
            system(['1dcat ' tsv_file_name '''[trans_x..rot_z_derivative1_power2]''' ' | tail -n 375 > ' subj_run '_regressors.1d']); %selects the last 375 lines of fields trans_x through rot_z...

            % Concatenate timepoints 8 through 382
            system(['3dTcat -prefix tcat_' subj_run '.nii ' file_name '''[8..382]''']); %8-382 corresponds to timepoints 9-383 in AFNI bc numbering starts at zero
            
            % Perform motion regression w/ AFNI
            system(['3dDeconvolve -input tcat_' subj_run '.nii ',...
                '-ortvec ' subj_run '_regressors.1d motion_' run ' ',...
                '-polort 3 -float -num_stimts 0 ',...
                '-fout -tout -x1D X_' subj_run '.xmat.1D ',...
                '-fitts fitts_' file_name ' ',...
                '-errts r_' file_name ' -bucket stats_' file_name(1:end-4) ' ',...
                '-mask ' mask_file]);
            
           % Perform despiking
           system(['3dDespike -NEW -nomask -prefix d' reg_files{kk} ' ' reg_files{kk}]); % add nomask and ssave options

           % Perform smoothing
           system(['3dmerge -1blur_fwhm 6 -doall -prefix s' dr_files{kk} ' ' dr_files{kk}]);

           %
           
        end
    end; toc %end of subject loop
    
%     %% Change OpenMP settings to not hog up so many cpus in parallel
%     setenv('OMP_NUM_THREADS', '4'); %Set # of openmp threads for parallel
% 
%     %% Despike in parallel with AFNI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tic
%     disp('Performing despiking.')
%     %check which of the reg files need to be despiked
% %     for ll = 1:size(dr_files, 1) 
% %         if isfile(dr_files(ll,1))
% %             disp([subj ' run-' num2str(runs_to_use(ll)) ' already exists. Ommitting from despiking.'])
% %             reg_files(ll) = []; %deletes the row from reg_files if the file already exists           
% %         end
% %     end
%     
%     parfor ii = 1:size(reg_files, 1)
%         tic; system(['3dDespike -nomask -prefix d' reg_files{ii} ' ' reg_files{ii}]); toc% add nomask and ssave options
%         tic; system(['3dDespike -NEW -nomask -prefix d' reg_files{ii} ' ' reg_files{ii}]); toc% add nomask and ssave options
% 
%         %system(['3dDespike -NEW -q -nomask -ssave spike_' reg_files{ii}  ' -prefix d' reg_files{ii} ' ' reg_files{ii}]);
%     end
%     toc
%     
%      %% Smoothing in parallel with AFNI, 6mm FWHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tic
%     disp('Performing smoothing.')
%     %check which of the dr files need to be smoothed
% %     for ll = 1:size(sdr_files, 1) 
% %         if isfile(sdr_files(ll,1))
% %             disp([subj ' run-' num2str(runs_to_use(ll)) ' already exists. Ommitting from smoothing.'])
% %             dr_files(ll) = []; %deletes the row from dr_files if the file already exists
% %         end
% %     end
% 
%     parfor ii = 1:size(dr_files, 1)
%         system(['3dmerge -1blur_fwhm 6 -doall -prefix s' dr_files{ii} ' ' dr_files{ii}]);
%         tic; system(['3dBlurToFWHM -quiet -FWHM 6 -prefix s' dr_files{ii} ' -mask ' mask_files{ii} ' -input ' dr_files{ii}]); toc; %alternative smoothing program
%     end
%     toc
    
    %% Intermediary file clean up temp files
%     s_files = dir('sdr*.nii');
%     if size(s_files, 1) == 4
%         disp('All files preprocessed. Deleting intermediary files.')
%         delete dr_*.nii;
%         delete r_*.nii;
%         delete tcat_*.nii;
%     end
%     
    %% Move back to the root directory
    cd(parent_dir_path)
    
end %end to the status == 0 conditional statement
toc