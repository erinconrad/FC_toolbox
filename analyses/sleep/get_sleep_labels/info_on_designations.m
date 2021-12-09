function info_on_designations

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
erin_des_folder = [results_folder,'analysis/sleep/erin_designations/'];
spikes_folder = [results_folder,'all_out/'];
%data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Loop through all patients in the erin_designations folder
listing = dir([erin_des_folder,'*.mat']);

for l = 1:length(listing)
    
    
end


end