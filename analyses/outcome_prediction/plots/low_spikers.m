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

locs = out.all_soz_locs;
lats = out.all_soz_lats;
npts = length(locs);
names = out.all_names;
good_spikes = nan(npts,1);
nspikes = nan(npts,1);


for ip = 1:npts
    curr_name = names{ip};
    
    % find corresponding row of table
    r = strcmp(T.name,curr_name);
    
    assert(sum(r) == 1)
    
    % get nspikes
    nspikes(ip) = T.spikes_min_car_(r);
    
    % good spikes
    good_spikes(ip) = T.PPV_car_(r) > 0.7;
end

T = table(names,nspikes,good_spikes,locs,lats);

%ignore_pts = T.good_spikes == 0;
%T(ignore_pts,:) = [];




end