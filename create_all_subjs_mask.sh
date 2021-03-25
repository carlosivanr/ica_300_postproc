# Creates a union mask of all input masks as long as they all have the same dimensions

# This will not work because some masks have a different volume dimensions
3dmask_tool -input /export/research/analysis/human/jhouck/cobre06_65007/carlos_work/abcd_flx/06_data/derivatives/fmriprep/sub-NDARINV????????/ses-baselineYear1Arm1/func/union_mask.nii -prefix all_subjs_mask.nii -union

# Alternatively, if the union_mask.nii is not yet made, can make an all subjects masks with fmriprep output masks
3dmask_tool -input /export/research/analysis/human/jhouck/cobre06_65007/carlos_work/abcd_flx/06_data/derivatives/fmriprep/sub-NDARINV0B7UGM1D/ses-baselineYear1Arm1/func/sub-NDARINV????????_ses-baselineYear1Arm1_task-MID_run-?_space-MNI152NLin6Asym_desc-brain_mask.nii.gz -prefix all_subjs_mask -union
