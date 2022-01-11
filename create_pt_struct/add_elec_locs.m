function add_elec_locs

overwrite = 1;

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
            continue
        end
        
    end
    
    out_native = return_elec_locs(pt_name,elec_path);
    out_mni = return_mni_locs(pt_name,elec_path);
    
    for i = 1:length(out_mni)
        out_mni(i).anatomy = out_native(i).anatomy;
    end
    
    pt(p).elecs = out_mni;
    save([data_folder,'pt.mat'],'pt');
end


end