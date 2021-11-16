function run_sleep_analyses

%% Parameters
main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
main_lats = {'Left','Right'};
main_soz = {'SOZ','Not SOZ'};
main{1} = main_locs;
main{2} = main_lats;
main{3} = main_soz;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);


end