function add_elec_locs

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
ieeg_folder = locations.ieeg_folder;
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
    out = return_elec_locs(pt_name,elec_path);
    pt(p).elecs = out;
    save([data_folder,'pt.mat'],'pt');
end


end