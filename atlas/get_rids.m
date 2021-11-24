function get_rids

%% Parameters
% which atlas

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
atlas_folder = [locations.main_folder,'data/atlas/'];
%parcellation_folder = [atlas_folder,'atlas_localizations/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load file containing conversion
T = readtable([atlas_folder,'HUPID_to_RID.csv']);

for p = 1:length(pt)
    
    %% Get name
    name = pt(p).name;
    %if strcmp(name,'HUP099'), name = 'HUP99'; end

    %% Convert HUPID to RID (this is the naming used for atlases)
    table_hup_ids = T.ieegportalsubjno;
    table_rids = T.record_id;
    match = contains(table_hup_ids,name);
    
    if sum(match) ~= 1, error('what'); end
    
    rid = table_rids(match);
    pt(p).rid = rid;
    %row = find(match);


    

end

save([data_folder,'pt.mat'],'pt');

end