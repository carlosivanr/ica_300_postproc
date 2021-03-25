%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carlos Rodriguez, Ph.D. Mind Research Network
% 03/22/2021
% Fix sub-NDARINVT2JJ42KE
% 
% This script was created to fix run-1 of the restin state data
% from ABCD study subject sub-NDARINVT2JJ42KE. This particular run
% only had 382 timepoints, whereas all other runs had 383 time points.
% Rather than excluding this subject, the last time point was 
% estimated as the mean of the last 10 TRs, for each voxel.

% This script is desinged to run in the subjects session 1 func
% directory containing the original zipped fmriprep file.

% Dependencies. Relies on AFNI 3dinfo to quickly count TRs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% -------------------Fix the time point-----------------------------
% Unzip run-1. Only run-1 needs to be interpolated
gunzip('sub-NDARINVT2JJ42KE_ses-baselineYear1Arm1_task-rest_run-1_space-MNIPediatricAsym_cohort-3_res-2_desc-preproc_bold_brain.nii.gz')
file = dir('sub-NDARINV*.nii');

% check time points
system(['3dinfo -n4 ' file.name])

% load the volume with SPM12 
volname = file.name; %file name of the rs-fmri run
v = spm_vol(volname); %spm command to create a variable with the header information
data = spm_read_vols(v); %spm command to read the entire volume with the corresponding header
%data = data(:,:,:,9:end); %excludes the first 8 volumes ie, dummy scans from Siemens ABCD data
[x,y,z,t] = size(data); %grabs the size dimensions of the data matrix x,y,z and time
tp = (1:t)'; %number of time points in a column vector, used for plotting during script development

% create the insert
ins = zeros(x,y,z,1); %initialize empty matrix for residuals
for xx = 1:size(data,1) % X
    for yy = 1:size(data,2) % Y
        for zz = 1:size(data,3) % Z axial slices
            Y = data(xx,yy,zz,:); %takes each voxel at all time points
            Y = Y(:); %converts Y into a column vector
            Y = Y(end-10:end);
            r = mean(Y); %performs the regression may need to use [a, R2, residual] = icatb_regress(y, X) script in the ICA Tool Box
            ins(xx,yy,zz,1) = r; %writes the residuals into the matrix
        end
    end
end

% write 3d volumes
V = v(1); %v is the header information from spm_vol function on the original image file
V.descrip = V(1).descrip; %modify description
for l = 1:size(data,4) %
    outputname = ['temp_' volname(1:end-4) '_' sprintf('%03d', tp(l)) '.nii']; %concatenates temp prefix, file name, image number, and .nii as output name
    V.fname = outputname; %sets up the file names with a prefix and an image number to take the 4D data to a series of 3D images
    spm_write_vol(V, data(:,:,:,l)); % Writes the new images
end

% write the last timepoint
V = v(1); %v is the header information from spm_vol function on the original image file
V.descrip = V(1).descrip; %modify description
outputname = ['temp_' volname(1:end-4) '_' num2str(383) '.nii']; %concatenates temp prefix, file name, image number, and .nii as output name
V.fname = outputname; %sets up the file names with a prefix and an image number to take the 4D data to a series of 3D images
spm_write_vol(V, ins(:,:,:)); % Writes the new images

    
% Merge 3d files back to one 4d file, this will overwrite the unzipped
% run-1 data, but is intended to be this way since the ica_300 will not
% unzip the original 382 time points data in this configuration.
nii_in = dir(['temp_' file.name(1:end-4) '_*.nii']);
nii_file = {};

for bb = 1:length(nii_in)
    nii_file{bb} = fullfile(pwd,nii_in(bb).name);
end

spm_file_merge(nii_file,volname,0,.8);

% Check time points
system(['3dinfo -n4 ' file.name])

% Clean up the temp files
delete(['temp_' file.name(1:end-4) '_*.nii']); %

% copy the dataset with the interpolated timepoint to the ica_input
% directory, so the original does not get unzipped.
copyfile(file.name, '/export/research/analysis/human/jhouck/abcd/ica_300/ica_input/sub-NDARINVT2JJ42KE');


%% -------------------Fix the .tsv file -----------------------------
% The tsv file will need an additional row of data to accomodate the
% added time point.

% Get the name of the .tsv file
tsv_file = dir('sub-*run-1*_regressors.tsv');

% conditional statement to prevent using an over-written .tsv to be renamed
% as original.
if isfile(['original_' tsv_file.name])
    disp('Original file already renamed and copied.')
else   
    copyfile(tsv_file.name, ['original_' tsv_file.name]);    
end

% read in the original .tsv file
temp_tsv = readtable(['original_' tsv_file.name], 'Filetype', 'text'); 

% modify the first 135 columns, with mean
for ii = 1:135
    if isa(temp_tsv.(ii), 'double') == 0
        value = mean(str2double(temp_tsv{373:382, ii})); %convert string to double, to calculate mean
        temp_tsv{383, ii} = {num2str(value)}; %replace with mean as string
     else   
        temp_tsv(383, ii) = {mean(temp_tsv{373:382, ii})}; %replace with the mean 
    end
end

% The step above fills in all numeric columns with zeros, so motion outliers and steady state cols
% don't need any modification

% 144 - 167
for ii = 144:167
    if isa(temp_tsv.(ii), 'double')
        temp_tsv(383, ii) = {mean(temp_tsv{373:382, ii})}; %replace with the mean as double
    else
        value = mean(str2double(temp_tsv{373:382, ii})); %convert string to double, to calculate mean
        temp_tsv{383, ii} = {num2str(value)}; %replace with mean as string
    end

end


% modify columns 136 through 143, aroma
for ii = 215:238
    temp_tsv(383, ii) = {mean(temp_tsv{end-10:end, ii})}; %replace with a zero
end

% Write to TSV
%Matlab currently has no way of writing .tsv files directly, but writing to
%.txt with \t delimiter and renaming to .tsv extension is a working hack
writetable(temp_tsv, [tsv_file.name(1:end-4) '.txt'], 'Delimiter', '\t');
movefile([tsv_file.name(1:end-4) '.txt'], [tsv_file.name(1:end-4) '.tsv']); %This overwrites the fmriprep output file

% The resulting file will have problems with wc -l, because it will return
% all lines, 384 including header, but the original output of fmriprep .tsv
% files will output 383 an not count the header.