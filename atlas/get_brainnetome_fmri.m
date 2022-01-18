function get_brainnetome_fmri_and_sc

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


%% load SC
d = dir([brainnetome_folder,'sc/BNA_SC_3D_246/*.nii.gz']);
SC = nan(length(d));
for j = 1:length(d)
    fprintf('SC ROI %d\n',j)
    nii = niftiread(fullfile(d(j).folder,d(j).name));
    roi_conn = PARCELLATE_IMAGE(nii,atlas_img); % get sc for that ROI
    SC(j,:) = roi_conn;
end
SC = (SC + SC')/2; % symmetrize
SC(logical(eye(length(d)))) = nan;

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
FC(logical(eye(length(d)))) = nan;

%% Save
out.FC = FC;
out.SC = SC;
save([out_folder,'brainnetome_sc_fc.mat'],'out');


end