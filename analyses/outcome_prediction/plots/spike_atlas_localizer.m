function spike_atlas_localizer

locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/data/'];

%% Load out file and get roc stuff
out = load([out_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
rl = out.all_rl;
soz = out.all_soz_bin;
loc = out.all_soz_locs;
npts = length(soz);

%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);




end