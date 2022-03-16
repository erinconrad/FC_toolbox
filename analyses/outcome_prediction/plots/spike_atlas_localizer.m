function spike_atlas_localizer

%% Parameters
which_atlas = 'aal_bernabei';

%% File locations
locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
spike_folder = [results_folder,'analysis/outcome/data/'];
atlas_folder = [results_folder,'analysis/atlas/'];


%% Load out file and get roc stuff
out = load([spike_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
rl = out.all_rl;
soz = out.all_soz_bin;
loc = out.all_soz_locs;
npts = length(soz);

%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);

%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

%% get atlas stuff
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;

end