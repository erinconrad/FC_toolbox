function surgical_info

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
script_folder = locations.script_folder;

%% Get pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Addpath
addpath(genpath(script_folder));

name = {};
surg = {};
resection_lat = {};
resection_loc = {};
ablation_lat = {};
ablation_loc = {};

for ip = 1:length(pt)
    name = [name;pt(ip).name];
    surg = [surg;pt(ip).clinical.surgery];
    resection_lat = [resection_lat;pt(ip).clinical.resection_lat];
    resection_loc = [resection_loc;pt(ip).clinical.resection_loc];
    ablation_lat = [ablation_lat;pt(ip).clinical.ablation_lat];
    ablation_loc = [ablation_loc;pt(ip).clinical.ablation_loc];
    
end

T = table(name,surg,resection_lat,resection_loc,ablation_lat,ablation_loc)

end