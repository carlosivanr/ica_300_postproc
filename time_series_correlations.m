function [subj_dir_path] = time_series_correlations(subj_dir_path)

%% Correlate Two Pipelines
% Carlos Rodriguez, Ph.D., Mind Research Network
% This function was written to correlated the time series of two seperate
% post fmrip prep pipelines: 1) composed of strictly AFNI commands to
% perform file truncation, regression, detrend, despike, and smoothing with
% openmp, 2) icatb regress, detrend, old afni 3dDepsike, and SPM
% smoothing.

%% Usage from matlab command line:
% argument in = absolute path to the subject directory
% ex.
% time_series_correlations('/export/reasearch/analysis/human/jhouck/abcd/ica_300/sub-NDAR*')

%
%% Usage from bash terminal:
% matlab -nodisplay -nojvm -r "icatb_ica_300_postproc('/absolute/path/to/subject-specific-fmriprep-output-directory');exit;"
%
%% Verify Directories, Files, etc.
ses = 'ses-baselineYear1Arm1';
[parent_dir_path, subj, ~] = fileparts(subj_dir_path); %retrieves parent directory and subject ID


%check for the directory existence
if isfolder([subj_dir_path filesep ses filesep 'func'])
    disp('Checking for output files to correlate.')
    cd([subj_dir_path filesep ses filesep 'func'])
    %list afni pipeline sdr files
    sdr = dir('sdr_sub*.nii');
    
    %list icatb pipeline sdr files
    isdr = dir('sdr_icatb*.nii');

    if size(sdr,1) == size(isdr,1)
        disp('Number of sdr and sdr_icatb files match, proceeding to check files dimensions.')
        for kk = 1:size(sdr, 1)
            % load the volume with SPM12
            % load the sdr (afni preprocessed data)
                s_volname = sdr(kk).name; %file name of the smoothed rs-fmri run
                sv = spm_vol(s_volname); %spm command to create a variable with the header information
                s_data = spm_read_vols(sv); %spm command to read the entire volume with the corresponding header

             % load the icatb (icatb/old despike data) using the sdr file
             % name to ensure runs are matched when correlating
                i_volname = ['sdr_icatb_' s_volname(5:end)]; %file name of the rs-fmri run
                iv = spm_vol(i_volname); %spm command to create a variable with the header information
                i_data = spm_read_vols(iv); %spm command to read the entire volume with the corresponding header

            % check file dimensions before calculating correlations
            if size(i_data) == size(s_data)
                disp('File dimensions match, proceeding to calculating correlations')
         
                % set the dimensions
                    [x,y,z,t] = size(i_data); %grabs the size dimensions of the data matrix x,y,z and time

                % calculate correlations
                    correlations = zeros(x,y,z); %initialize empty matrix for residuals
                    for xx = 1:size(i_data,1) % X
                        for yy = 1:size(i_data,2) % Y
                            for zz = 1:size(i_data,3) % Z axial slices
                                Y = i_data(xx,yy,zz,:); %takes each voxel at all time points for the icatb data
                                Y = Y(:);
                                X = s_data(xx,yy,zz,:); %takes each voxel at all time points for the afni data
                                X = X(:);
                                r = corrcoef(X,Y);
                                correlations(xx,yy,zz) = r(2,1); %writes the residuals into the matrix
                            end
                        end
                    end
                    
                 % Calculate the average correlation
                 
                 
                 % write out the correlations volume
                    disp('Writing to a correlation .nii file.')
                    V = iv(1); %v is the header information from spm_vol function on the original image file
                    V.descrip = [V(1).descrip ' Correlations between AFNI and ICATB post fmri prep processing pipelines']; %modify description
                    outputname = ['corr_' i_volname(11:end)]; %concatenates temp prefix, file name, image number, and .nii as output name
                    V.fname = outputname; %sets up the file names with a prefix and an image number to take the 4D data to a series of 3D images
                    spm_write_vol(V, correlations(:,:,:)); % Writes the new images 
             else
                disp('File dimensions do no match. Check file output.')
                return
            end
        end 
    else
        disp('Number of sdr and sdr_icatb files do not match. Check output.')
        return
    end

else
    disp('/func directory not found. Check file path.')
    return
end
