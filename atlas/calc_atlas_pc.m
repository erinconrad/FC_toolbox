function avg_pc = calc_atlas_pc

%% To do
%{
- remove soz or resected ch
- campbell wrote script in matlab to go from MNI to AAL
%}

%% Parameters
fs = 200;

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
atlas_folder = [data_folder,'atlas/'];
atlas_file = [atlas_folder,'HUP_atlas_final.mat'];

%% Load the atlas file
load(atlas_file);

avg_pc = pc_vector_calc(eeg,fs,2);
avg_pc = wrap_or_unwrap_adjacency_fc_toolbox(avg_pc);

end