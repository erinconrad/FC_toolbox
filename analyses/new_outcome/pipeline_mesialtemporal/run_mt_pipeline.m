%% RUN MT pipeline
function run_mt_pipeline(whichPts)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_path = [results_folder,'edf_out/'];
data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load pt folder
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

for i = 1:length(whichPts)
    ip = whichPts(i);
    name = pt(ip).name;
    mt_patient_stitch(edf_path,name)

end


end