function add_rids_to_struct

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load table with RIDs
T = readtable('Manual validation.xlsx','Sheet','RIDs');

for p = 1:length(pt)
    
    %% Get name
    name = pt(p).name;

    %% Find row in table with this name
    r = strcmp(T.name,name);
    assert(sum(r)==1) % make sure you can find it

    %% Get the corresponding RID
    rid = T.RIDs(r);

    %% Add it to the pt struct
    pt(p).rid = rid;

    %% Save pt file
    save([data_folder,'pt.mat'],'pt');


end

end