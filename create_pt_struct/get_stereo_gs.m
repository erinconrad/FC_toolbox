function get_stereo_gs

%% Get file locs
locations = fc_toolbox_locs;
data_folder = [locations.main_folder,'data/'];
script_folder = locations.script_folder;

%% Get pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Addpath
addpath(genpath(script_folder));

%% Get the xls file with the file start times
T = readtable('Manual validation.xlsx','Sheet','Stereo');

%% Loop over patients
for p = 1:length(pt)
    name = pt(p).name;
    
    % find corresponding row in table
    row = strcmp(name,T.name);
    assert(sum(row) == 1);
    
    stereo_str = T.G_SOrStereo{row};
    if strcmp(stereo_str,'Stereo')
        pt(p).clinical.stereo = 1;
    elseif strcmp(stereo_str,'G/S')
        pt(p).clinical.stereo = 0;
    else
        error('what')
    end
    
end


save([data_folder,'pt.mat'],'pt');

end