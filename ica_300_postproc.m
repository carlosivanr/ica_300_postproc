%% ICA 300 post fmriprep processing function
% Carlos Rodriguez, Ph.D. 
% This function was created for the ica_300 project to unzip the best of 4 
% runs of resting state fmri from 300 participants of the ABCD study.
% Voxelwise timecourses are detrended using gift icatb_detren
% Motion parameters are then regressed out of the voxelwise timecourses
% using gift icatb_regress. 
%
% Dependencies include Gift ICA toolbox, SPM12, and AFNI set up to work
% through bash. Matlab parallel pool for smoothing multiple subjects at 
% once is also a dependency. Run this script from the root directory 
% containing all of the subjects fmriprep output.
%
% 12/28/20 - Updated 3dTcat, parallel processing for smoothing
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
%
%
% Modifications to make:
% - convert to function, using the subject directory as input
% - modify script to not require a subject loop
% - modify script to perform QA on it's own with out qa_merge
% - after w/ in subject QA, could try a batch process to perform detrend,
% regression, despike, and smooth

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change directory to ica_300 fmriprep output directory
%cd('/export/research/analysis/human/jhouck/abcd/ica_300/fmriprep_out/fmriprep');%uncomment when ready to deploy on ica_300

%***% 
        rdr = '/export/research/analysis/human/jhouck/cobre06_65007/carlos_work/ica_300/ica_300_concatenate_runs'; %sets the root directory, during development

%% Load the QA merged data to figure out which files to unzip in cases of more than 4 resting state runs
% The line below is to be uncommented when finished to apply to the entire data set
%runs_to_use = readtable('/export/research/analysis/human/jhouck/cobre06_65007/carlos_work/ica_300/ica_300/ica_300_fmriprep_qa.csv');

%****% Development of script relies on a truncated runs-to-use .csv file
load('dev_runs_to_use.mat') % To be commented when deploying to full ica_300.

%% Initialize batch for parallel smoothing with SPM
%  for smoothing parfor loops batch needs to be a 1 x n cell array. Each 
%  cell should be a 1 x 2 cell containing the matlabbatch data
%batch = cell(1, size(runs_to_use,1)*4);
batch = {}; %initialize empty cell 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through each subject in the qa merge file

for ii = 2:size(runs_to_use, 1)
    % Initialize variables within each subject loop iteration
    % Initialize subject
    subj = runs_to_use.bids_id{ii};
    
    % Change into the subject directory
    cd([subj '/ses-baselineYear1Arm1/func']); % uses relative path, not absolute
    
    % List all of rest fmri zipped files
    files = dir('sub-*_desc-preproc_bold.nii.gz'); %lists the rest fmri files

    % Use the qa merge file to determine which runs to unzip if more than 4
    qc_files = runs_to_use.runs_to_use_if_more_than_4{ii};

    % Unzipping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If the qc files is empty that means, then there are only 4 runs to
    % unzip
    disp('Checking files to unzip.')
    if isempty(qc_files) 
        for jj = 1:size(files,1)
            filezp = files(jj).name; % temporary string of zipped file names to access the string
            fileun = files(jj).name(1:end-3); %unzipped file with the string stripped of the .gz for if then statement
            run = files(jj).name(strfind(files(jj).name, 'run'):strfind(files(jj).name, 'run')+4);
            if isfile(fileun)
                disp([subj ' ' run ' file unzipped.'])
            else
                disp([subj ' ' run ' file is zipped. Unzipping...'])
                gunzip(filezp);
            end            
        end 
    % If the qc files is not empty, then there were more than 4 runs and
    % the best 4 runs need to be unzipped according to qa merge .csv file 
    else 
        ind = sscanf(qc_files, '%g,', [4, inf]).' ; % converts string to number array to be used as an index of the string of files to unzip
        for jj = 1:size(files,1) %this loop will use ind to index the best 4 files, based off of the least number of motion outliers from fmriprep
            filezp = files(ind(jj)).name; %see above, although here, ind serves as th index
            fileun = files(ind(jj)).name(1:end-3); %see above
            run = files(jj).name(strfind(files(jj).name, 'run'):strfind(files(jj).name, 'run')+4); %identifies which run

            if isfile(fileun)
                disp([subj ' ' run ' file unzipped.'])
            else
                disp([subj ' ' run ' file is zipped. Unzipping...'])
                gunzip(filezp);
            end
        end        
    end
    
    % Detrend, regress motion, and despike %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % list the unzipped resting state bold files
    disp('Performing detrend and regression of voxelwise timecourses.')
    rs_files = dir('sub*bold.nii');
    
    % List confounds.tsv files in subject's directory
    cons = dir('*-confounds*.tsv'); 

    % for each rs fmri file that has been unzipped
    for kk = 1:size(rs_files, 1) % 1 through size rs_files along the 1st dim
        run = rs_files(kk).name(strfind(rs_files(kk).name, 'run'):strfind(rs_files(kk).name, 'run')+4); %identifies which run
        
%         Older code using tdfread. Takes a bit longer than using readtable.        
%         tic
%         % load corresponding confounds_regressors.tsv file
%         temp_tsv = tdfread(cons(kk).name, '\t'); %n.b. Some data types will be in char, results in struct
%         
%         % Prep the tsv file, replace with 1st cell nan with zeros, select only the columns
%         % needed for motion regression
%         
%         %Create index of fields that need conversion of char to num using OR operators
%         names = fieldnames(temp_tsv); %list out the field names
%         
%         %prep the derivative columns because they have nans in the first
%         %field and are in character instead of num format
%         cols = contains(names,'derivative'); %Index of all of the derivative measures
%         col_names = names(cols); %Cell that contains indexed field names
%             
%         %Convert char to num and replace initial nan cells with 0
%         for ll = 1:size(col_names,1)
%             temp_tsv.(col_names{ll})(1,:) = '0'; %[]; %Delete the first cell that contains the n/a of derivatives
%             temp_tsv.(col_names{ll}) = str2num(temp_tsv.(col_names{ll})); %Converts char to num. n.b. str2double did not work properly
%             %temp_tsv.(col_names{ll}) = [0; temp_tsv.(col_names{ll})]; % concatenates zero with the rest of the values %[mean(c.(names_ind{i})(1:end)); c.(names_ind{i})]; %Concatenates the average and remaining data
%         end
%         toc
        
        %% Insertion Point for readtable of .tsv files
        tic
        
        % load corresponding confounds_regressors.tsv file
        temp_tsv = readtable(cons(kk).name, 'Filetype', 'text');  % reads the .tsv files, but still results in processing
        
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
        toc
        
        %% 
        % Gather the prepped motion parameters for motion regression
        cols = contains(names,'trans') | contains(names, 'rot'); % index of trans and rot motion parameters
        col_names = names(cols); %Cell that contains indexed field names
    
        %Convert struct to table for preparing regressor design matrix X
        %regressors = struct2table(temp_tsv); %Convert structure to table
        clear temp_tsv %clear temp_tsv from memory
        
        % Regression design matrix
        rX = temp_tsv(:, col_names);
        %rX = regressors(:, col_names); %selects the column names from regressors, used in older temp_tsv code section
        rX = rX{:,:}; % selects only the values and not variable names for regression
        
        % load the volume with SPM12
        volname = rs_files(kk).name; %file name of the rs-fmri run
        v = spm_vol(volname); %spm command to create a variable with the header information
        data = spm_read_vols(v); %spm command to read the entire volume with the corresponding header
        [x,y,z,t] = size(data); %grabs the size dimensions of the data matrix x,y,z and time
        
        % perform motion regression on each voxel's time course
        tp = (1:size(data, 4))'; %number of time points in a column vector, used for plotting during script development
        
        % Using icatb functions, takes about 17 minutes for one file
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

        % Write the residuals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
               
        % Merge 3d files to one 4d file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Could test fslmerge, but will leave files gzipped. Will need to
            % unzip or figure out how to change the fsl environmental variable 
            % from NIFT_GZ to NIFTI
            % command = ['fslmerge -t r_' volname ' temp_*.nii -tr .8'];
            % system(command)
            % delete temp_*.nii;
         
        disp(['Merging 3d files for ' run '.'])
        nii_in = dir('temp_*.nii');
        nii_file = {};
        for bb = 1:length(nii_in)
            nii_file{bb} = fullfile(pwd,nii_in(bb).name);
        end
        spm_file_merge(nii_file,['r_' volname],0,.8);

        % Clean up the temp files
        delete temp_*.nii;
    
        % Despike %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each file that has been motion regressed
        % perform the 3d Despike
        disp(['Performing despiking for ' run '.'])
%%%        %command = ['3dDespike -prefix dr_' volname ' -nomask -ssave
        %spike_' volname ' r_' volname]; % to be deployed when fully
        %test
        command = ['3dDespike -prefix dr_' volname ' -nomask r_' volname]; %doesn't save spike data while testing

        system(command)
        delete r_*.nii;

%         % Convert AFNI .brik & .head to .nii %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         disp(['Converting to NIFTI for ' run '.'])
%         command = ['3dAFNItoNIFTI -prefix dr_' volname ' d+orig.BRIK'];
%         system(command)
%         delete d+*;

    end %ends the loop for each resting state fmri run
        
    % Set up parallel for loop for smoothing with SPM12 and 6mm FWHM %%%%%%
    disp('Setting up parallel loop batch for smoothing.')
    dr_files = dir('dr*.nii'); %list the despiked and regressed files
    for kk = 1:size(dr_files, 1)
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'smooth_files';
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {
                                                                     {
                                                                     [rdr '/' subj '/ses-baselineYear1Arm1/func/' dr_files(kk).name]
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
    
    % Smoothing in parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Proceeding to smoothing.')
    parfor xx = 1:size(batch,2)
         try
                out{xx} = spm_jobman('run',batch{xx})
         catch
                 out{xx} = 'failed';
          end
    end
   
    % Move back to the root directory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(rdr)

end %ends the loop for each subject

disp('Post-processing finished in ...')
toc
