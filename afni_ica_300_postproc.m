function [subj_dir_path] = afni_ica_300_postproc(subj_dir_path)

%% ICA 300 post fmriprep processing function
% This function was created for the ica_300 project to unzip the best of 4
% runs of resting state fmri from 300 participants of the ABCD study.
% Motion parameters are regressed out of the voxelwise timecourses,
% residual time courses are then despiked and smoothed.

% Dependencies include: AFNI version 20.2.02
% Set AFNI_AUTOGZIP = NO in .afnirc.

% Areas for further development;


%% Usage from matlab command line:
% argument in = absolute path to the subject directory
% ex.
% afni_ica_300_postproc('/export/reasearch/analysis/human/jhouck/abcd/ica_300/sub-NDAR*')

%% Usage from bash terminal:
% matlab -nodisplay -nojvm -r "afni_ica_300_postproc('absolute/path/to/subject-specific-fmriprep-output-directory');exit;"

%% Set log features
tic
[parent_dir_path, subj, ~] = fileparts(subj_dir_path); %retrieves parent directory and subject ID
format_out = 'mm-dd-yyyy';
proc_date = datestr(now, format_out);

if exist('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/logs', 'dir') == 0
    mkdir('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/logs')
end
    
diary(['/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/logs/' subj '_' proc_date '_log.txt'])
disp(datestr(now, 'mm-dd-yyyy HH:MM:SS PM'))
output_dir = '/export/research/analysis/human/jhouck/abcd/ica_300/ica_input';

%% Set environmental variables and check AFNI Version
setenv('PATH', [getenv('PATH') ':/export/research/analysis/human/jhouck/shared/tools/abin_20202']); %Set AFNI Version by appending to the path
setenv('OMP_NUM_THREADS', '16'); %Set # of openmp threads

% Check AFNI Version
[~, ver] = system('afni -ver');
disp(ver)
ver_num = str2double(ver(58+5:58+6)); %numeric of the first two digits of the AFNI version number

if ver_num >= 20
    disp('AFNI version set correctly set.')
else
    disp('Error: AFNI Version 20.2.02 or higher required. Exiting..')
    return
end

%% Verify directories, files, and errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
disp(['Processing ' subj])
disp(subj_dir_path)
fprintf('\n')
ses = 'ses-baselineYear1Arm1';
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

            %status of grep for "No errors found!", 0=no error/successful, 1=error may be reported, 2=permission denied for grep
            if contains(cmdout, 'No errors') %status == 0
                disp('No errors reported in .html file. Proceeding to post fmriprep processing.')
            elseif status == 1 %grep for "No errors" returns no results, check for errors
                [~, cmdout] = system(['grep "Errors" ' parent_dir_path filesep subj '.html']); % search the html report
                if contains(cmdout, 'Errors') %
                    disp('Errors in fmriprep html report found. Exiting ica_300_postproc.')
                    log_msg = [subj ', Errors in fmriprep html report reported'];
                    fid = fopen('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/post_proc_log_file.txt', 'a'); %open the log file to append message
                    fprintf(fid, '%s, %s\n', datestr(now, 0), log_msg);
                    fclose(fid);
                    return
                end
            elseif status == 2
                disp('func directory exists and at least 4 rs-fMRI files found, but fmriprep html report is unverifiable. Possible permission error.')
                log_msg = [subj ', fmriprep html report unverifiable. Possible permission error'];
                fid = fopen('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/post_proc_log_file.txt', 'a'); %open the log file to append message
                fprintf(fid, '%s, %s\n', datestr(now, 0), log_msg);
                fclose(fid);
                return
            end
        else
            disp('No fmriprep html report found, errors unverifiable. Exiting ica_300_postproc.')
                log_msg = [subj ', fmriprep html report not found'];
                fid = fopen('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/post_proc_log_file.txt', 'a'); %open the log file to append message
                fprintf(fid, '%s, %s\n', datestr(now, 0), log_msg);
                fclose(fid);
                return
        end
    else
        %disp([ num2str(size(rs_runs,1)) ' resting state files found. Exiting ica_300_postproc.'])
        fprintf(2, [ num2str(size(rs_runs,1)) ' resting state files found. Exiting ica_300_postproc.\n']) %prints error message in red
        log_msg = [subj ', '  num2str(size(rs_runs,1)) ' resting state files found.'];
        fid = fopen('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/post_proc_log_file.txt', 'a');
        fprintf(fid, '%s, %s\n', datestr(now, 0), log_msg); %open the log file to append message
        fclose(fid);
        return

    end
else
    disp('No func directory found. Exiting ica_300_postproc.')
    log_msg = [subj ', No func directory found.'  ];
    fid = fopen('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/post_proc_log_file.txt', 'a');
    fprintf(fid, '%s, %s\n', datestr(now, 0), log_msg); %open the log file to append message
    fclose(fid);
    return
end
clear status

%% Conditional statement to proceed if func directory exists, at least 4 rs-fMRI runs exist, and no error reported by fmriprep
if contains(cmdout, 'No errors to report!')
    
    % Check if subject output directory exists, if not, create on
    if exist([output_dir filesep subj], 'dir') == 0
        mkdir([output_dir filesep subj])       
    end
    
    % Change directory to ica_300 fmriprep subj directory
    cd([subj_dir_path filesep ses filesep 'func']);%
    

    % QC number of runs, timepoints, and number of motion outliers
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(rs_runs,1) >= 4
        fprintf('\n')
        disp('Checking which files to process.')
        tsv_files = dir('sub*_regressors.tsv'); %list .tsv confound regressor files
        num_outliers = zeros(1, size(tsv_files, 1)); %preallocate number of motion outliers per run
        num_timepoints = zeros(1, size(tsv_files, 1)); %preallocate number of time points per run
        for ii = 1:size(tsv_files,1)
            [~, mot] = system(['grep -o "motion_outlier" ' tsv_files(ii).name ' | wc -l']);
            num_outliers(ii) = str2double(mot);
            %[~, ntp] = system(['wc -l <' tsv_files(ii).name]); %This count method does not include header on orig .tsv files
            [~, ntp] = system(['cut -c 1 ' tsv_files(ii).name ' | wc -l']); %This count method includes header
            %num_timepoints(ii) = str2double(ntp);
            num_timepoints(ii) = str2double(ntp)-1; %subtracts one from bc of the header
        end
        run_info = [1:size(tsv_files,1); num_timepoints; num_outliers]'; % Create a matrix with the run_info
        T = array2table(run_info);  % Write run info to a .csv file as part of logging information regarding file processing
        T.Properties.VariableNames = {'run', 'timepoints', 'motion_outliers'};
        writetable(T, [output_dir filesep subj filesep subj '_run_info.csv']); %writes out a table with run diagnostics
        disp(T); clear T
        run_info(run_info(:, 2) ~= 383, :) = [] ; %deletes any rows with timepoints not equal to 383.
    end

    % If remaining runs with 383 timepoints are less than 4, exit out of the function and write to error log
    if size(run_info, 1) < 4
        disp('Not enough runs to process. Number of timepoints may be less than expected. Exiting ica_300_postproc.')
        log_msg = [subj ', Not of enough runs to process. Number of timepoints may be less than expected.'];
        fid = fopen('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/post_proc_log_file.txt', 'a');
        fprintf(fid, '%s, %s\n', datestr(now, 0), log_msg); %open the log file to append message
        fclose(fid);
        return
    elseif size(run_info, 1) >= 4 %rank the remaining files according to the number of motion outliers
        run_info = sortrows(run_info, 3); %sort run_info according to the number of motion outliers
        runs_to_use = run_info(1:4, 1); %select the the first 4 runs which contain the least number of motion outliers
        runs_to_use = sort(runs_to_use); %sort runs_to_use to reflect run acquisition
        disp(['Using runs ' num2str(runs_to_use(1)) ', ' num2str(runs_to_use(2)) ', ' num2str(runs_to_use(3)) ', and ' num2str(runs_to_use(4)) '.'])
        fprintf('\n')
    end

    % Copy and unzip resting state files
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:4 %this loop uses runs_to_use as an index to unzip files, because they may not always be in order
        filezp = rs_runs(runs_to_use(jj)).name; %name of the file when zipped
        fileun = filezp(1:end-3); %name of the file when unzipped  to check if it exists.
        run = rs_runs(runs_to_use(jj)).name(strfind(rs_runs(runs_to_use(jj)).name, 'run'):strfind(rs_runs(runs_to_use(jj)).name, 'run')+4); %identifies which run
        if isfile([output_dir filesep subj filesep fileun]) %does the file exist in the output directory?
            disp([run ' rest file unzipped.'])
        else
           %disp([run ' rest file not found in subject output directory. Copying and unzipping...'])
           %system(['cp ' filezp ' ' output_dir '/' subj '/']);
           %system(['gunzip ' output_dir '/' subj '/' filezp]);
           
           disp([run ' rest file not found in subject output directory. Unzipping to output directory...'])
           system(['gunzip -k -c ' filezp ' > ' output_dir '/' subj '/' fileun]);

        end
    end
    fprintf('\n')
    
    % Copy and unzip mask files
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:4 %this loop uses runs_to_use as an index to unzip files
        filezp = masks(runs_to_use(jj)).name; %name of the file when zipped
        fileun = filezp(1:end-3); %name of the file when unzipped  to check if it exists.
        run = ['run-' num2str(runs_to_use(jj))];
        if isfile([output_dir filesep subj filesep fileun])
            disp([run ' mask file unzipped.'])
        else
           %disp([run ' mask file not found in subject output directory. Copying and unzipping...'])
           %system(['cp ' filezp ' ' output_dir '/' subj '/']);
           %system(['gunzip ' output_dir '/' subj '/' filezp]);
           
           disp([run ' rest file not found in subject output directory. Unzipping to output directory...'])
           system(['gunzip -k -c ' filezp ' > ' output_dir '/' subj '/' fileun]);
        end
    end
    fprintf('\n')
    
    % Copy .tsv files
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Copying .tsv files.')
    for jj = 1:4 %this loop uses runs_to_use as an index to unzip files
        system(['cp ' subj '_' ses '_task-rest_run-' num2str(runs_to_use(jj)) '_desc-confounds_regressors.tsv ' output_dir '/' subj '/']);
    end
    fprintf('\n')   
    
    clear run_info num_timepoints num_outliers ntp mot fileun filezp masks

    % Post fmriprep processing: file truncation, motion regression, despiking, and smoothing
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd([output_dir '/' subj]);
    % Initialize empty cells for file checking
    reg_files = cell(4,1);
    dr_files = cell(4,1);
    sdr_files = cell(4,1);
    mask_files = cell(4,1);

    %tic;
    for kk = 1:4
        run = ['run-' num2str(runs_to_use(kk))]; %identifies which run is processing
        file_name = [subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %rs-fmri file to concatenate
        mask_file = [subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-brain_mask.nii'];%mask file name
        tsv_file_name = [subj '_' ses '_task-rest_' run '_desc-confounds_regressors.tsv']; %concatenates a .tsv file name

        reg_files{kk,1} = ['r_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for despiking regressed files
        dr_files{kk,1} = ['dr_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for smoothing despiked files
        sdr_files{kk,1} = ['sdr_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for smoothing despiked files
        mask_files{kk,1} = mask_file; %initialize a cell for despiking

        % check if the final smoothed file exists
        if isfile(sdr_files{kk}) ~= 1
            %disp([run ' final output exists. Skipping...'])
            %return
        %else
            %fprintf('\n') %This line could be deleted 
            disp('Performing post fmriprep processing.')
            fprintf('\n')
            disp([subj ' ' run '. Truncating files, extracting regressors, performing motion regression, and despiking.'])
            fprintf('\n')
            % Set up a string variable for file output
            subj_run = [subj '_' run];

            % Prepare a 1d file with the motion regressors 24 columns
            % representing the 6 motion parameters, derivatives, squares of
            % parameters, and squares of derivatives
            system(['1dcat ' tsv_file_name '''[trans_x..rot_z_derivative1_power2]''' ' | tail -n 375 > ' subj_run '_regressors.1D']); %selects the last 375 lines of fields trans_x through rot_z...

            % Some confounds files are arranged differently
            [~, cons] = system(['awk ' '''{print NF}'' ' subj_run '_regressors.1D | sort -nu | tail -n 1']); %counts the number of columns
            cons = str2double(cons);
            if cons < 24
                %fprintf('\n')
                disp('Number of regressors is less than 24. Trying a different configuration of extracting regressors.')
                delete([subj_run '_regressors.1D']);
                %delete([subj_run '_regressors.1D']);
                system(['1dcat ' tsv_file_name '''[trans_x..rot_z_power2]''' ' | tail -n 375 > ' subj_run '_regressors.1D']); %selects the last 375 lines of fields trans_x through rot_z...
                [~, cons] = system(['awk ' '''{print NF}'' ' subj_run '_regressors.1D | sort -nu | tail -n 1']); %counts the number of columns, multiple apostrophes used to escape '{print NF}' and add a space afterward
                cons = str2double(cons);
                fprintf('\n')
                if cons ~= 24
                    fprintf('\n')
                    disp(['Error in extracting motion regressors for ' subj ' consider reviewing confounds.tsv file and modify code. A total of 24 motion regressors were not extracted.'])
                    fprintf('\n')
                    log_msg = [subj ', Error in extracting motion regressors. Review confounds .tsv file.'];
                    fid = fopen('/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/post_proc_log_file.txt', 'a');
                    fprintf(fid, '%s, %s\n', datestr(now, 0), log_msg); %open the log file to append message
                    fclose(fid);                                        
                    return
                end

            end

            % Select and Demean motion parameters
            % columns of motion parameters [1:4:24] -1; %columns 1 through 24 by 4s, minus 1 bc afni indexing starts at 0
            system(['1dcat -sel ''[0,4,8,12,16,20]'' ' subj_run '_regressors.1D > ' subj_run '_motion_regressors.1D']); %selects the motion parameters only

            % Demeans the translation and rotation motion_regressors
            system(['1d_tool.py -infile ' subj_run '_motion_regressors.1D -set_nruns 1 -demean -write ' subj_run '_demeaned_regressors.1D']); %set to demean the motion parameters

            % Derivatives & Squares
            system(['1dcat -sel ''[1,2,3,5,6,7,9,10,11,13,14,15,17,18,19,21,22,23]'' ' subj_run '_regressors.1D > ' subj_run '_motion_derivatives.1D']); %

            % Demeans the derivatives
            system(['1d_tool.py -infile ' subj_run '_motion_derivatives.1D -set_nruns 1 -demean -write ' subj_run '_demeaned_derivatives.1D']); %set to demean the motion parameters

            % Concatenate timepoints 8 through 382
            system(['3dTcat -prefix tcat_' subj_run '.nii ' file_name '''[8..382]''']); %AFNI numbering starts at zero, therefore 8:382 corresponds to timepoints 9:383 bc they're all siemens files

            % Perform motion regression w/ AFNI 3dDeconvolve
            system(['3dDeconvolve -input tcat_' subj_run '.nii ',...
                '-ortvec ' subj_run '_demeaned_regressors.1D motion_' run ' ',...
                '-ortvec ' subj_run '_demeaned_derivatives.1D motion_' run ' ',...
                '-polort 3 -float -num_stimts 0 ',...
                '-fout -tout -x1D X_' subj_run '.xmat.1D ',...
                '-xjpeg X_' subj_run '.jpg ',...
                '-fitts fitts_' file_name ' ',...
                '-errts r_' file_name ' -bucket stats_' file_name(1:end-4) ' ',...
                '-mask ' mask_file]);

           % Perform despiking
           disp('Performing despiking.')
           system(['3dDespike -NEW -nomask -prefix d' reg_files{kk} ' ' reg_files{kk}]); % add nomask and ssave options

           % Perform smoothing
           disp('Performing smoothing.')
           system(['3dmerge -1blur_fwhm 6 -doall -prefix s' dr_files{kk} ' ' dr_files{kk}]);
        else
            disp(['sdr files for ' run ' were already processed and skipped.'])

        end
    end %toc %end of subject loop

    %% Intermediary file clean up temp files
    fprintf('\n')
    s_files = dir('sdr*.nii');
    if size(s_files, 1) == 4
        disp('All files processed. Deleting intermediary files.')
        delete dr_*.nii;
        delete r_*.nii;
        delete tcat_*.nii;
    else
        disp('Not all runs processed. Check error log.')
    end

    %% Move back to the root directory
    cd(parent_dir_path)

end %end to the status == 0 conditional statement
time = toc;
if time < 60
    disp(['Elapsed processing time for ' subj ' was ' num2str(time) ' seconds.'])
else
    disp(['Elapsed processing time for ' subj ' was ' num2str(time/60) ' minutes.'])
end
diary off
