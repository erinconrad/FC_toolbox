function get_brainnetome_fmri

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/atlas/'];
data_folder = [locations.main_folder,'data/'];
brainnetome_folder = [data_folder,'atlas/brainnetome/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% load atlas
atlas_img = niftiread([brainnetome_folder,'BN_Atlas_246_2mm.nii.gz']);


%% load FC
% Get list of regions
d = dir([brainnetome_folder,'fc/BNA_FC_3D_246/*.nii.gz']);
FC = nan(length(d));
% Loop over regions
for j = 1:length(d)
    fprintf('FC ROI %d\n',j)
    nii = niftiread(fullfile(d(j).folder,d(j).name));
    roi_conn = PARCELLATE_IMAGE(nii,atlas_img); % get sc for that ROI
    FC(j,:) = roi_conn;
end
FC = (FC + FC')/2; % symmetrize
FC(~eye(length(d))) = nan;

%% Save
out.FC = FC;
save([out_folder,'brainnetome_FC.mat'],'out');


end