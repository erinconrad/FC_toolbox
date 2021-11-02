function count_missing_sz_pts

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
data_folder = [locations.main_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

missing_idx = [];
missing_names = {};

portal_idx = [];
portal_names = {};

for p = 1:length(pt)
    if ~isfield(pt(p).ieeg.file(1),'sz_times')
        missing_idx = [missing_idx;p];
        missing_names = [missing_names;pt(p).name];
        continue
    end
    
    if strcmp(pt(p).ieeg.file(1).sz_time_source,'portal times, worst source')
        portal_idx = [portal_idx;p];
        portal_names = [portal_names;pt(p).name];
    end
    
end

fprintf('\nMissing seizure times:\n');
table(missing_idx,missing_names)

fprintf('\nPortal seizure times:\n');
table(portal_idx,portal_names)

end