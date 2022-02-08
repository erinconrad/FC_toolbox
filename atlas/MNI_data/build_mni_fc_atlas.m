function build_mni_fc_atlas

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/atlas/'];
int_folder = [results_folder,'analysis/intermediate/'];
data_folder = [locations.main_folder,'data/'];
mni_folder = [data_folder,'MNI_open_ieeg/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load mni data
mni = load([mni_folder,'MatlabFile.mat']);

%% Get the eeg data
eeg = mni.Data_W;

end