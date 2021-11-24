function build_fc_atlas

%% Parameters
net_m = 2;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
data_folder = [locations.main_folder,'data/'];
spikes_folder = [results_folder,'all_out/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load pt folder
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

names = {};

%% Loop over patients
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    name = summ.name;
    
    %% Find corresponding pt index
    found_it = 0;
    for ip = 1:length(pt)
        if strcmp(pt(ip).name,name)
            found_it = 1;
            break
        end
    end
    if ~found_it, error('what'); end
    
    %% get rid and labels
    rid = pt(ip).rid;
    elabels = summ.labels;
    
    %% Get atlas
    out = get_atlas_parcellations(rid,elabels,name);
    
    %% Skip if empty
    if isempty(out.enum)
        continue
    end
    
    %% Load the spike file
    fname = [spikes_folder,name,'_pc.mat'];
    if ~exist(fname,'file')
        fprintf('\nCannot find spike file for %s, skipping...\n',name);
        continue
    end
    
    %% reconcile files (deal with changes in electrode names)
    out = net_over_time(pc,pt,ip);
    out = reconcile_files(out);
    
    %% Get fc
    fc = out.montage(net_m).net;
    fc = wrap_or_unwrap_adjacency_fc_toolbox(fc);
    
    %% Average over time
    fc = nanmean(fc,3);
    
    
end

end