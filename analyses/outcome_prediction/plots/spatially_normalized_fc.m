function spatially_normalized_fc


locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
atlas_folder = [results_folder,'analysis/atlas/'];

%% Load out file
out = load([out_folder,'main_out.mat']);
out = out.out;

%% get stuff
soz = out.all_soz_bin;
loc = out.all_soz_locs;
fc = out.all_fc;
npts = length(soz);
labels = out.all_labels;

fit_distance_model(all_locs,all_conn)

end