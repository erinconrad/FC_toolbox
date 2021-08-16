function add_elec_locs

overwrite = 0;

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
script_folder = locations.script_folder;

%% Get pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Addpath
addpath(genpath(script_folder));
box_path = locations.box_folder;
elec_path = [box_path,'CNT Implant Reconstructions/'];

%% Loop over patients and get locs
for p = 1:length(pt)
    pt_name = pt(p).name;
    
    if isfield(pt(p),'elecs') == 1 && ...
            ~isempty(pt(p).elecs)
        if overwrite == 0
            fprintf('\nskipping %s\n',pt(p).name);
        end
        continue
    end
    
    out = return_elec_locs(pt_name,elec_path);
    pt(p).elecs = out;
    save([data_folder,'pt.mat'],'pt');
end


end