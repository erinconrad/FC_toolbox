%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [locations.main_folder,'data/'];
edf_path = [results_folder,'edf_summ_out/'];
sleep_stage_path = [results_folder,'edf_out/'];
out_folder = [results_folder,'analysis/new_outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Load pt folder to get names
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
npts = length(pt);
all_names = cell(npts,1);
for ip = 1:npts
    all_names{ip} = pt(ip).name;
end

%% Loop over patients
for p = 1:npts

    name = all_names{p};


    % Load the edf summary file
    if exist([edf_path,name,'/summ.mat'],'file') == 0, continue; end
    info = load([edf_path,name,'/summ.mat']);
    info = info.out;
    times = info.all_times;


    % also load the sleep stages
    stage = load([sleep_stage_path,name,'/sleep_stage.mat']);
    sout = stage.sout;


    % Get times of sleep transitions
    st_dates= sout.Summary(2:end,2);
    st_times = sout.Summary(2:end,3);
    stage = sout.Summary(2:end,4);

    % convert to seconds in ieeg
    seeg_secs = convert_transitions_to_ieeg_secs(st_dates,st_times);

    % get the sleep stages for appropriate times
    [~,seeg_stage] = match_seeg_times(times(:,1),ones(size(times,1),1),seeg_secs,stage);

end