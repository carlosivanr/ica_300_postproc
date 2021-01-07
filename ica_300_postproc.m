function [subj_dir_path] = ica_300_postproc(subj_dir_path)

%% ICA 300 post fmriprep processing function
% Carlos Rodriguez, Ph.D. 
% This function was created for the ica_300 project to unzip the best of 4 
% runs of resting state fmri from 300 participants of the ABCD study.
% Voxelwise timecourses are detrended using gift icatb_detren
% Motion parameters are then regressed out of the voxelwise timecourses
% using gift icatb_regress. 
%
% Dependencies include Gift ICA toolbox, SPM12, and AFNI set up to work
% through bash. Set AFNI_AUTOGZIP = NO in .afnirc. Set OpenMP threads,
% e.g 'export OMP_NUM_THREADS=16'. Set FSL to NIFTI output, e.g.
% export FSLOUTPUTTYPE=NIFTI. Matlab parallel pool for 
% smoothing multiple subjects at once is also a dependency.
%
%% Usage
% argument in = absolute path to the subject directory
% ex.
% ica_300_postproc('/export/reasearch/analysis/human/jhouck/abcd/ica_300/sub-NDAR*')
% requires write permissions for error_logging in the directory from which
% it is called.
%
% To run from bash:
% matlab -nodisplay -nojvm -r "ica_300_postproc('absolute/path/to/subject-specific-fmriprep-output-directory');exit;"
%
%% Revision History
% 12/28/20 - Updated 3dTcat and parallel processing for smoothing
% 12/29/20 - Removed 3dTcat for concatenating, added code to read .tsv 
%   confound regressors. Added detrend and regress functions, file merge of
%   3d volumes, 3dDespike
% 12/30/20 - Added displays, cleaned up file merge code section. Modified
% code to run from subjects 2 through 8 in testing.
%   check timepoints after regression
%   check slices in 3dDespike
% 01/04/21 - Modified reading of .tsv files to use readtable instead of
% tdfread which speeds up the process, but still requires conversion of
% 'n/a's to double. Set AFNI 3dDespike prefix to include .nii to bypass the
% .brik to .nii conversion step. Added -nomask and -ssave options. 
% Smoothing set to run within each subject folder.
% 01/05/21 - Modified script to work as a matlab that takes the subject 
% specific directory as an argument. Changed preproc_bold to
% preproc_brain_bold.nii files to process which are masked
% 01/06/21 - Added parfor loop functionality to 3dDespike, added
% functionality to use fslmerge if FSLOUTPUTTYPE is set to NIFTI.
%
% Modifications to make:
% - convert to function, using the subject directory as input
% - modify script to not require a subject loop
% - modify script to perform QA on it's own with out qa_merge
% - after w/ in subject QA, could try a batch process to perform detrend,
% regression, despike, and smooth
% - add session input argument, in case there will be session 2 data to
% process
% - parfor loop for unzipping files
%
% Testing to perform:
% - run function from matlab
% - run function from bash
%
% Loose ends to tie up
% - ensure regress options are detrend and regress
% - ensure -nomask and -ssave are reset after testing
% - ensure delete temp* and r_* files are reset

%% Code section for testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subjects = dir('sub-NDAR*'); %list all of the files and folders that begin with sub-* prefix
% subjects = subjects([subjects.isdir]); %modify list to contain only directories
% subj_dir_path  = [subjects(1).folder filesep subjects(1).name]; %set subj_dir_path to subject 01

%% Verify directories, files, and errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ses = 'ses-baselineYear1Arm1';
[parent_dir_path, subj, ~] = fileparts(subj_dir_path); %retrieves parent directory and subject ID
        % [~, parent_folder_name] = fileparts(parent_dir_path) ;
        % parent_folder_name = 'Dir'

disp('Checking func directory ...')
% check if the func directory exists
if isfolder([subj_dir_path filesep ses filesep 'func'])
    disp('func directory found. Checking files...')
    
    % check the number of rs files
    rs_runs = dir([subj_dir_path filesep ses '/func/*preproc_bold_brain.nii.gz']); %list all the zipped rest runs
    if size(rs_runs, 1) >= 4
        disp('At least 4 resting state files found. Checking fmriprep html report for errors...')
        
        % check the html report for errors
        if isfile([parent_dir_path '/' subj '.html']) %if file exists, proceed to the next lines
            %T.exists_html(i) = 1; %fill in the table with the .html exists
            [status, cmdout] = system(['grep "No errors to report!" ' parent_dir_path filesep subj '.html']); % check the existence of the html report
            
            %status of grep for "No errors found!", 0=no error/successful, 1=not found, 2=permission denied for grep            
            if contains(cmdout, 'No errors') %status == 0
                disp('No errors found. Proceeding to postprocessing.')
            elseif status == 1
                disp('fmriprep html report not found.')
                log_msg = [subj ', fmriprep html report not found'];
                fid = fopen('post_proc_log_file.txt', 'a'); %open the log file to append message
                fprintf(fid, '%s: %s\n', datestr(now, 0), log_msg);
                fclose(fid);
            elseif status == 2
                disp('func directory exists and at least 4 rs-fMRI files found, but fmriprep html report is unverifiable due to permission denied.')
                log_msg = [subj ', fmriprep html report unverifiable, possible permission error'];
                fid = fopen('post_proc_log_file.txt', 'a'); %open the log file to append message
                fprintf(fid, '%s: %s\n', datestr(now, 0), log_msg);
                fclose(fid);
            end
        else
            disp('No fmriprep html report found, errors unverifiable. Exiting ica_300_postproc')
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
end

%% Conditional statement to proceed if func directory exists, at least 4 rs-fMRI runs exist, and no error reported by fmriprep
    if contains(cmdout, 'No errors')
    % Change directory to ica_300 fmriprep output directory
    cd([subj_dir_path filesep ses filesep 'func']);%

    % Determine which runs should be used once in the sub/ses/func director. %%%%%%%%%%%%%%%
    % If more than 4 resting state runs exist, use the runs with the least
    % number of motion outliers from the .tsv files
    if size(rs_runs,1) > 4
        tsv_files = dir('*_regressors.tsv'); %list .tsv confound regressor files
        n_m_outliers = zeros(1, size(tsv_files,1)); %preallocate number of motion outliers per run
        for ii = 1:size(tsv_files,1)
            [~, cmdout] = system(['grep -o "motion_outlier" ' tsv_files(ii).name ' | wc -l']);
            n_m_outliers(ii) = str2double(cmdout);
        end
        n_m_outliers(2,:) = 1:size(tsv_files,1); %add an index for the rs_run
        n_m_outliers = sortrows(n_m_outliers'); %sort the rows descending and transpose
        runs_to_use = n_m_outliers(1:4, 2); %Select the runs with the least amount of motion outliers
    else
        runs_to_use = '';
    end

    % Unzip Files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If runs_to_use is empty, then there are only 4 runs to unzip
    disp('Checking files to unzip.')
    if isempty(runs_to_use) %if runs to use is empty, then only 4 runs of rs-fmri exist
        for jj = 1:size(rs_runs,1)
            filezp = rs_runs(jj).name; %temporary string of zipped file names to access the string
            fileun = filezp(1:end-3); %unzipped file with the string stripped of the .gz for if then statement
            run = rs_runs(jj).name(strfind(rs_runs(jj).name, 'run'):strfind(rs_runs(jj).name, 'run')+4);
            if isfile(fileun)
                disp([run ' file unzipped.'])
            else
                disp([run ' file is zipped. Unzipping...']) 
                system(['gunzip -k ' filezp]);
            end
        end
        runs_to_use = 1:4; %set runs to use for indexing tsv files

    % If runs_to_use is not empty, then there were more than 4 runs and
        % the best 4 runs need to be unzipped 
        else 
            %ind = n_m_outliers(:,2); % 
            for jj = 1:4 %this loop will use n_m_outliers to index the best 4 files
                filezp = rs_runs(runs_to_use(jj)).name; %see above, although here, ind serves as the index
                fileun = filezp(1:end-3); %see above
                run = rs_runs(runs_to_use(jj)).name(strfind(rs_runs(runs_to_use(jj)).name, 'run'):strfind(rs_runs(runs_to_use(jj)).name, 'run')+4); %identifies which run
                if isfile(fileun)
                    disp([run ' file unzipped.'])
                else
                   disp([run ' file is zipped. Unzipping...']) 
                   system(['gunzip -k ' filezp]);
                end
            end
    end

    % Detrend, regress motion, and despike %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Performing detrend and regression of voxelwise timecourses.')
    uz_files = dir('sub*bold_brain.nii'); % list the unzipped resting state bold files

    % List confounds.tsv files in subject's directory
    tsv_files = dir('*_regressors.tsv'); %list .tsv confound regressor files
    tsv_files = tsv_files(sort(runs_to_use, 'ascend')); % sort the runs to use

    % for each rs fmri file that has been unzipped
    for kk = 1:size(uz_files, 1) % 1 through size uz_files along the 1st dim
        run = uz_files(kk).name(strfind(uz_files(kk).name, 'run'):strfind(uz_files(kk).name, 'run')+4); %identifies which run is processing
        subj = uz_files(kk).name(1:19); %select just the subject name
        tsv_file_name = [subj '_ses-baselineYear1Arm1_task-rest_' run '_desc-confounds_regressors.tsv']; 

        %load corresponding confounds_regressors.tsv file    
        temp_tsv = readtable(tsv_file_name, 'Filetype', 'text');

        % Prep the tsv file, replace 1st cell 'n/a's with zeros, select only the columns
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
        rX = rX{:,:}; % selects only the values and not variable names for regression
        clear temp_tsv %clear temp_tsv from memory

        % load the volume with SPM12
        volname = uz_files(kk).name; %file name of the rs-fmri run
        v = spm_vol(volname); %spm command to create a variable with the header information
        data = spm_read_vols(v); %spm command to read the entire volume with the corresponding header
        [x,y,z,t] = size(data); %grabs the size dimensions of the data matrix x,y,z and time

        % perform motion regression on each voxel's time course
        tp = (1:size(data, 4))'; %time points in a column vector, used for plotting and numbering 3d volumes when writing out .nii files
            %figure    
            %plot(tp, data(x,y,z,:))

        % Using icatb functions, takes about 17 minutes for one file
        residuals = zeros(x,y,z,t); %initialize empty matrix for residuals
        for xx = 1:size(data,1) % X
            for yy = 1:size(data,2) % Y
                for zz = 1:size(data,3) % Z axial slices
                    Y = data(xx,yy,zz,:); %takes each voxel at all time points
                    Y = Y(:); %converts Y into a column vector
%                    dY = icatb_detrend(Y, 1, size(data,4), 1); % Gift tool box function detrends linearly similar to matlab's native detrend function
%                    [~,~,r] = icatb_regress(dY, rX); %performs the regression may need to use [a, R2, residual] = icatb_regress(y, X) script in the ICA Tool Box
                    [~, ~, r] = icatb_regress(Y, rX); %performs the regression may need to use [a, R2, residual] = icatb_regress(y, X) script in the ICA Tool Box
                    residuals(xx,yy,zz,:) = r; %writes the residuals into the matrix

                end
            end
        end

        % Write the residuals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % excludes the first 8 volumes ie, dummy scans, for Siemens ABCD
        % study data only
        disp(['Writing residuals to 3d volumes for ' run '.'])
        V = v(1); %v is the header information from spm_vol function on the original image file
        V.descrip = [V(1).descrip ' Voxelwise time courses detrended and regressed with icatb_detrend and icatb_regress (GroupICATv4.0c)']; %modify description
        for l = 9:size(data,4) %starts at 9 to exclude dummy scans
            outputname = ['temp_' volname(1:end-4) '_' sprintf('%03d', tp(l)) '.nii']; %concatenates temp prefix, file name, image number, and .nii as output name
            V.fname = outputname; %sets up the file names with a prefix and an image number to take the 4D data to a series of 3D images
            spm_write_vol(V, residuals(:,:,:,l)); % Writes the new images
        end

        % Merge 3d files to one 4d file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['Merging 3d files for ' run '.'])
           
        % If FSL set to NIFT, use fslmerge, else use SPM to merge files    
        [~, cmdout] = system('echo $FSLOUTPUTTYPE'); % Check FSL output type
        if contains(cmdout, 'NIFTI')
            disp('Merging with FSL')
            system(['fslmerge -tr r_' volname ' temp_' volname(1:end-4) '*.nii' .8]); %uses volname with a r_ as a prefix to write out the file    
        else
            disp('Merging with SPM')
            nii_in = dir(['temp_' volname(1:end-4) '*.nii']); %lists temp prefix concatenated with volname to ensure not all temp* get merged
            nii_file = {};
            for bb = 1:length(nii_in)
                nii_file{bb} = fullfile(pwd,nii_in(bb).name);
            end
            spm_file_merge(nii_file,['r_' volname],0,.8);
        end
        
    end %ends the resting state file loop
    
    %% Despike %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each file that has been motion regressed
    % perform the 3d Despike
    %disp(['Performing despiking for ' run '.'])
    reg_files = dir('r_sub*bold_brain.nii');
    parfor ii = 1:size(reg_files, 1)
        system(['3dDespike -q -prefix d' reg_files(ii).name ' ' reg_files(ii).name]); % add nomask and ssave options
    end

    %% Set up parallel for loop for smoothing with SPM12 and 6mm FWHM %%%%%%%%%%
    disp('Setting up parallel loop batch for smoothing.')
    dr_files = dir('dr*brain.nii'); %list the despiked and regressed files
    batch = {}; %initialize empty cell 

    for kk = 1:size(dr_files, 1)
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'smooth_files';
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {
                                                                     {
                                                                     [parent_dir_path '/' subj '/ses-baselineYear1Arm1/func/' dr_files(kk).name]
                                                                     }
                                                                     }';
        matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Named File Selector: smooth_files(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{2}.spm.spatial.smooth.fwhm = [6 6 6]; %set FWHM smoothing kernel to 6mm
        matlabbatch{2}.spm.spatial.smooth.dtype = 0;
        matlabbatch{2}.spm.spatial.smooth.im = 0;
        matlabbatch{2}.spm.spatial.smooth.prefix = 's';
        subj_batch{1} = matlabbatch;
        batch = [batch subj_batch]; %appends to a growing list of batches
    end

    disp('Proceeding to smoothing.')
    parfor xx = 1:size(batch,2)
         try
                out{xx} = spm_jobman('run',batch{xx})
         catch
                 out{xx} = 'failed';
          end
    end

    %% Intermediary file clean up despiked and regressed files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delete dr_*.nii;
    delete r_*.nii;
    delete temp_*.nii:

    % Move back to the root directory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(parent_dir_path)
end %end to the status == 0 conditional statement