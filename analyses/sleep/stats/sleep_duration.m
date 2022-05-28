function sleep_duration(durations)

%% Locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
time_roc_folder = [out_folder,'time_roc/'];
if ~exist(time_roc_folder,'dir')
    mkdir(time_roc_folder)
end

end