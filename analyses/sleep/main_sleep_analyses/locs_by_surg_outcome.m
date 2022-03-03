function out = locs_by_surg_outcome(disc)

%% Parameters
which_locs = {'temporal','other cortex'};
which_lats = {'Left','Right'};

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/seizure_times/'];
int_folder = [results_folder,'analysis/backup_intermediate_Feb26_good_spikes/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% Prep array of spikes counts by location and sleep wake
spikes_loc = nan(npts,2,2); % temp/other cortex and then left/right
spikes_loc_ws = nan(npts,2,2,2); % wake-sleep
data_struct = {'left temporal','right temporal';'left other cortex','right other cortex'};

%% Loop over patients
for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get main things
    loc = summ.ana_loc;
    lat = summ.ana_lat;
    spikes = summ.spikes;
    ad = summ.ad;
    labels = summ.labels;
    
    %% Switch locs to temporal/other cortex
    loc = bin_into_temp_extra(loc);
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);

    loc = loc(~ekg,:);
    lat = lat(~ekg);
     
    spikes = spikes(~ekg,:); 
    labels = labels(~ekg);
     
    %% Determine "wake" and "sleep" times
    [sleep,wake] = find_sleep_wake(ad,[],disc);
     
    %% Bin spikes by loc
    for i = 1:2 % temp, other cortex
        for j = 1:2 % left, right
            
            idx = strcmp(loc,which_locs{i}) & strcmp(lat,which_lats{j});
            
            spikes_loc(p,i,j) = nanmean(spikes(idx,:),'all');
            spikes_loc_ws(p,i,j,1) = nanmean(spikes(idx,wake),'all');
            spikes_loc_ws(p,i,j,2) = nanmean(spikes(idx,sleep),'all');
            
        end
    end
    
end

out.spikes_loc = spikes_loc;
out.spikes_loc_ws = spikes_loc_ws;
out.data_struct = data_struct;

end