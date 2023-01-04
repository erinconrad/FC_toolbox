%% RUN MT pipeline
function run_mt_pipeline(whichPts,overwrite)

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_path = [results_folder,'edf_out/'];
edf_summ_path = [results_folder,'edf_summ_out/'];
data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load pt folder
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

if isempty(whichPts)
    whichPts = 1:length(pt);
end

for i = 1:length(whichPts)
    ip = whichPts(i);
    name = pt(ip).name;
    mt_patient_stitch(edf_path,edf_summ_path,name,overwrite)

end


end