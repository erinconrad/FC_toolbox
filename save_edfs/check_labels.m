function check_labels

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
out_dir = [results_folder,'edf_out/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Add sleep seeg folder to path
addpath(genpath(locations.sleep_seeg_folder))

% Loop over folders in edf out
listing = dir(out_dir);
for l = 1:length(listing)
    
    if exist([out_dir,listing(l).name,'/sleep_stage.mat']) == 0
        continue
    end

    

    sout = load([out_dir,listing(l).name,'/sleep_stage.mat']);
    sout = sout.sout;

    if ~isfield(sout,'name')
        continue
    end

    sout.name
    sout.ChL
    pause

end

end