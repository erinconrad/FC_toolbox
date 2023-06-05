%% RUN MT pipeline
function run_mt_pipeline(whichPts,overwrite)

%% Parameters
attempt_remove_oob = 0; % try to remove electrodes out of the brain

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_path = [results_folder,'edf_out/'];
edf_summ_path = [results_folder,'edf_summ_out/'];
data_folder = [locations.main_folder,'data/'];

if ~exist(edf_summ_path,'dir')
    mkdir(edf_summ_path)
end

%% Create an overlap log file
overlap_log_file = [edf_summ_path,'overlap_log_file.csv'];
if exist(overlap_log_file,'file') == 0
    T = table([],[],[],[],[],[],[],[],[],...
        'VariableNames',{'name','edf_file','sz_start_time','sz_end_times',...
        'sz_time_aligned','sz_end_time_aligned','attempted_start_time','attempted_end_time','num_attempts'});
    writetable(T,overlap_log_file);
end


% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% LOad validation file
validation_file = [scripts_folder,'spike_detector/Manual validation.xlsx'];
szT = readtable(validation_file,'Sheet','SOZ');
mT = readtable(validation_file,'Sheet','strange_elec_names');


%% Load pt folder
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

if isempty(whichPts)
    whichPts = 1:length(pt);
end

for i = 1:length(whichPts)
    ip = whichPts(i);
    name = pt(ip).name;
    mt_patient_stitch(pt,ip,edf_path,edf_summ_path,name,overwrite,overlap_log_file,szT,mT,attempt_remove_oob);

end


end