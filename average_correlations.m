% This script was created to calculate the average correlation and standard
% deviation of the corr .nii files which contain correlations between an
% older icatb based preprocessing pipeline and a newer afni based pipeline.
% Dependencies, SPM12

subjects = dir('sub-*');
subjects = subjects([subjects.isdir]); %modify subjects to contain only directories
output = [];

for ii = 1:size(subjects, 1)
    subj = subjects(ii).name;
    cd ([subj '/ses-baselineYear1Arm1/func/']);
    
    corr_files = dir('corr_*.nii');
    
    T = table('Size', [4, 7], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'subj', 'run', 'mean_r', 'std', 'median', 'min', 'max'});
    
    for jj = 1:size(corr_files,1)
        run = corr_files(jj).name(strfind(corr_files(jj).name, 'run-'):strfind(corr_files(jj).name, 'run-')+4); %identify where the position of run- string begins and adds 4 to extract which run is being read
        
        %Read the volume
        volname = corr_files(jj).name; %file name of the smoothed rs-fmri run
        v = spm_vol(volname); %spm command to create a variable with the header information
        data = spm_read_vols(v); %spm command to read the entire volume with the corresponding header
        
        % Prep the data and calculate summary stats
        data_vec = data(:);
        data_vec = rmmissing(data_vec); %remove missing entries, i.e. nans        
        
        % Fill in the table
        T.subj(jj) = subj;
        T.run(jj) = run;
        T.mean_r(jj) = mean(data_vec);
        T.std(jj) = std(data_vec);
        T.median(jj) = median(data_vec);
        T.min(jj) = min(data_vec);
        T.max(jj) = max(data_vec);
        
        %writetable(T, [subj '_corr_info.csv']); %writes out a table with
        %average correlations and standard deviations
       
    end
    disp(T)
   
    % Concatenate table contents
    output = [output; table2array(T)];
    
    cd ../../..
    
end

%Write a table with all of the output
output = array2table(output);
output.Properties.VariableNames = {'subj', 'run', 'mean_r', 'std', 'median', 'min', 'max'};
writetable(output, 'pipeline_corr_info.csv'); %writes out a table to csv
