function low_spikers

locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/data/'];

%% Load out file and get roc stuff
out = load([out_folder,'main_out.mat']);
out = out.out;

%% LOad manual validation file
T = readtable('manual validation.xlsx');


end