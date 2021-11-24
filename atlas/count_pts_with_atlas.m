function count_pts_with_atlas

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
data_folder = [locations.main_folder,'data/'];
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
    labels = pt(ip).labels;
    
    %% Get atlas
    out = get_atlas_parcellations(rid,elabels);

    if ~isempty(out.enum)
        names = [names;name];
    end
    
end

table(names)
fprintf('\nThere are %d of %d patients with some atlas data.\n',length(names),npts);


end