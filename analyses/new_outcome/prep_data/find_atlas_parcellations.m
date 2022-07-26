function out = find_atlas_parcellations(locs,atlas)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/atlas/'];
data_folder = [locations.main_folder,'data/'];
atlas_folder = [data_folder,'atlas/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

switch atlas
    case 'aal'
        [atlas_names,atlas_nums] = aal_region_to_name([atlas_folder,'AAL116_WM.txt'],[]);
        [~, mni_roi, ~] = nifti_values_erin(locs,[atlas_folder,'AAL116_WM.nii']);
        out.enum = mni_roi;
        out.atlas_names = atlas_names';
        out.atlas_nums = atlas_nums;
        
        
        
    case 'brainnetome'
        library = readtable([atlas_folder,atlas,'_library.xlsx']);
        atlas_nums = library.nums;
        atlas_names = library.names;
        mni_roi = mni2atlas(locs,[atlas_folder,'brainnetome/BN_Atlas_246_2mm.nii.gz']);
        out.enum = mni_roi;
        out.atlas_names = atlas_names;
        out.atlas_nums = atlas_nums;
end

%% Assign atlas names to electrodes
out.enames = cell(size(locs,1),1);
[lia,locb] = ismember(out.enum,out.atlas_nums);
out.enames(lia) = out.atlas_names(locb(lia));


end