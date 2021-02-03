function [subj_dir_path] = icatb_ica_300_postproc(subj_dir_path)

%% ICA 300 post fmriprep processing function based on icatb functions
% Carlos Rodriguez, Ph.D., Mind Research Nework
% Created to process fmriprep output of ABCD study fmri data

%% Usage from matlab command line:
% argument in = absolute path to the subject directory
% ex.
% icatb_ica_300_postproc('/export/reasearch/analysis/human/jhouck/abcd/ica_300/sub-NDAR*')
% requires write permissions for error_logging in the directory from which
% it is called.
%
%% Usage from bash terminal:
% matlab -nodisplay -nojvm -r "icatb_ica_300_postproc('/absolute/path/to/subject-specific-fmriprep-output-directory');exit;"
%
tic
%% Verify Directories, Files, etc.
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
        T.Properties.VariableNames = {'run', 'timepoints', 'motion_outliers'};
        writetable(T, [subj '_run_info.csv']); %writes out a table with run diagnostics
        disp(T); clear T
        run_info(run_info(:, 2) ~= 383, :) = [] ; %deletes any rows with timepoints not equal to 383.
    end
    
    % If remaining runs with 383 timepoints are less than 4, exit out of the function and write to error log
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
    
    % Unzip resting state files
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    clear run_info num_timepoints num_outliers ntp mot fileun filezp rs_runs masks
    
    % Post fmriprep processing: file truncation, motion regression, despiking, and smoothing
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Performing post fmriprep processing.')
    
    % Initialize empty cells for file checking
    %reg_files = cell(4,1);
    dr_files = cell(4,1);
    sdr_files = cell(4,1);
    %mask_files = cell(4,1);
    
    % for loop process the resting state runs
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk = 1:4
        run = ['run-' num2str(runs_to_use(kk))]; %identifies which run is processing
        file_name = [subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %rs-fmri file to concatenate
        %mask_file = [subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-brain_mask.nii'];%mask file name
        tsv_file_name = [subj '_' ses '_task-rest_' run '_desc-confounds_regressors.tsv']; %concatenates a .tsv file name
        
        reg_files{kk,1} = ['r_icatb_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for despiking regressed files
        dr_files{kk,1} = ['dr_icatb_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for smoothing despiked files
        sdr_files{kk,1} = ['sdr_icatb_' subj '_' ses '_task-rest_' run '_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii']; %initialize a cell for smoothing despiked files
        %mask_files{kk,1} = mask_file; %initialize a cell for despiking

        % check if the final smoothed file exists
        if isfile(sdr_files{kk})
            disp([run ' final output exists. Skipping...'])
        else
            disp([subj ' ' run ' final output does not exist. Extracting regressors and performing regression of motion parameters.'])
            
            % Set up a string variable for file output
            % subj_run = [subj '_' run];
            
            % Prepare the motion regressors
            % load corresponding confounds_regressors.tsv file
            temp_tsv = readtable(tsv_file_name, 'Filetype', 'text');  % reads the .tsv files, but still results in processing
        
            % Prep the tsv file, replace with 1st cell nan with zeros, select only the columns
            % needed for motion regression

            %Create index of fields that need conversion of char to num using OR operators
            names = fieldnames(temp_tsv); %list out the field names

            %prep the derivative columns because they have 'n/a's in the first
            %field and are in character instead of num format
            cols = contains(names,'derivative'); %Index of all of the derivative measures
            col_names = names(cols); %Cell that contains indexed field names

            %Replace initial 'n/a's with 0 and convert char to double and 
            for ll = 1:size(col_names,1)
                temp_tsv.(col_names{ll})(1,:) = {'0'}; %replace 'n/a's wit zero char
                temp_tsv.(col_names{ll}) = str2double(temp_tsv.(col_names{ll})); %Converts char to double
            end
 
            % Gather the prepped motion parameters for motion regression
            cols = contains(names,'trans') | contains(names, 'rot'); % index of trans and rot motion parameters
            col_names = names(cols); %Cell that contains indexed field names
         
            % Regression design matrix
            rX = temp_tsv(:, col_names);
            rX = rX{9:end,:}; % selects only the values and not variable names for regression
            
           % load the volume with SPM12 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            volname = file_name; %file name of the rs-fmri run
            v = spm_vol(volname); %spm command to create a variable with the header information
            data = spm_read_vols(v); %spm command to read the entire volume with the corresponding header
            data = data(:,:,:,9:end); %excludes the first 8 volumes ie, dummy scans from Siemens ABCD data
            [x,y,z,t] = size(data); %grabs the size dimensions of the data matrix x,y,z and time
            tp = (1:t)'; %number of time points in a column vector, used for plotting during script development
            
            % perform motion regression on each voxel's time course 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Using icatb functions, takes about 17 minutes for one file
            
            if isfile(reg_files{kk}) == 0 %if the icatb_file already exists, don't run the regression again
                residuals = zeros(x,y,z,t); %initialize empty matrix for residuals
                for xx = 1:size(data,1) % X
                    for yy = 1:size(data,2) % Y
                        for zz = 1:size(data,3) % Z axial slices
                            Y = data(xx,yy,zz,:); %takes each voxel at all time points
                            Y = Y(:); %converts Y into a column vector
                            dY = icatb_detrend(Y, 1, size(data,4), 1); % Gift tool box function detrends linearly similar to matlab's native detrend function
                            [~,~,r] = icatb_regress(dY, rX); %performs the regression may need to use [a, R2, residual] = icatb_regress(y, X) script in the ICA Tool Box
                            residuals(xx,yy,zz,:) = r; %writes the residuals into the matrix

                        end
                    end
                end
            end

            % Write the residuals 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['Writing residuals to 3d volumes for ' run '.'])
            V = v(1); %v is the header information from spm_vol function on the original image file
            V.descrip = [V(1).descrip ' Voxelwise time courses detrended and regressed with icatb_detrend and icatb_regress (GroupICATv4.0c)']; %modify description
            for l = 1:size(data,4) %
                outputname = ['temp_' volname(1:end-4) '_' sprintf('%03d', tp(l)) '.nii']; %concatenates temp prefix, file name, image number, and .nii as output name
                V.fname = outputname; %sets up the file names with a prefix and an image number to take the 4D data to a series of 3D images
                spm_write_vol(V, residuals(:,:,:,l)); % Writes the new images
            end

            % Merge 3d files back to one 4d file
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
            disp(['Merging 3d files for ' run '.'])
            nii_in = dir(['temp_' file_name(1:end-4) '_*.nii']);
            nii_file = {};
            
            for bb = 1:length(nii_in)
                nii_file{bb} = fullfile(pwd,nii_in(bb).name);
            end
            
            spm_file_merge(nii_file,['r_icatb_' volname],0,.8);

            % Despike with AFNI
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isfile(dr_files{kk}) %delete any pre existing file because AFNI does not overwrite
                delete(dr_files{kk});
            end
            
            disp(['Performing despiking for ' run '.'])            
            setenv('PATH', [getenv('PATH') ':/export/research/analysis/human/jhouck/shared/tools/abin_20202']);  %Set AFNI Version by appending to the path
            system(['3dDespike -prefix dr_icatb_' volname ' -nomask r_icatb_' volname])
         
            % Smooth wtih SPM
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isfile(dr_files{kk})
                clear matlabbatch         
                matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'despiked_files';
                matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{[subj_dir_path '/ses-baselineYear1Arm1/func/' dr_files{kk}]}};
                matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Named File Selector: despiked_files(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
                matlabbatch{2}.spm.spatial.smooth.fwhm = [6 6 6];
                matlabbatch{2}.spm.spatial.smooth.dtype = 0;
                matlabbatch{2}.spm.spatial.smooth.im = 0;
                matlabbatch{2}.spm.spatial.smooth.prefix = 's';       
                spm_jobman('run', matlabbatch); %process files in serial.
            end
            
            % Clean up the temp files
            delete(['temp_' file_name(1:end-4) '_*.nii']); %
            
        end %to the nested conditional statement checking for the sdr file

    end %ends the loop for each resting state fmri run
    
end % ends the main conditional checking for the status of file and directory checks

disp('Post-processing finished in ...')
toc
