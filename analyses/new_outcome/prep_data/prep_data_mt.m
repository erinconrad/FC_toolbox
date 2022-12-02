function prep_data_mt

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_path = [results_folder,'edf_out/'];
out_folder = [results_folder,'analysis/new_outcome/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the current output file
data = load([out_folder,'main_out.mat']);
data = data.out;
all_names = data.all_names;
npts = length(all_names);

% Goal is to build in same structure
% Structure is all vs wake vs sleep
all_coh = cell(npts,3);
all_bp = cell(npts,3);
all_spikes = cell(npts,3);
all_rl = cell(npts,3);
all_pearson = cell(npts,3);
all_labels = cell(npts,1);


%% Loop over patients
for p = 1:npts

    name = all_names{p};

    % Load the edf summary file
    if exist([edf_path,name,'/summ.mat'],'file') == 0, continue; end
    info = load([edf_path,name,'/summ.mat']);
    info = info.out;

    % get the time varying measurements
    spike_counts = info.all_spike_counts;
    rl = info.all_rl;
    bp = info.all_bp;
    coh = info.all_coh;
    pc = info.all_pc;
    times = info.all_times;
    labels = info.labels;
    all_labels{p} = labels;
    fs = info.fs;

    rl  = rl/fs;

    % also load the sleep stages
    stage = load([edf_path,name,'/sleep_stage.mat']);
    sout = stage.sout;

    % Get times of sleep transitions
    st_dates= sout.Summary(2:end,2);
    st_times = sout.Summary(2:end,3);
    stage = sout.Summary(2:end,4);

    % convert to seconds in ieeg
    seeg_secs = convert_transitions_to_ieeg_secs(st_dates,st_times);

    % get the sleep stages for appropriate times
    [~,seeg_stage] = match_seeg_times(times(:,1),ones(size(times,1),1),seeg_secs,stage);

    sleep = strcmp(seeg_stage,'N3') | strcmp(seeg_stage,'N2') | ...
        strcmp(seeg_stage,'N1') | strcmp(seeg_stage,'R');
    wake = strcmp(seeg_stage,'W');

    % Save all, wake, and sleep
    all_coh{p,1} = squeeze(nanmean(coh,1));
    all_coh{p,2} = squeeze(nanmean(coh(wake,:,:,:),1));
    all_coh{p,3} = squeeze(nanmean(coh(sleep,:,:,:),1));

    all_spikes{p,1} = squeeze(nanmean(spike_counts,1));
    all_spikes{p,2} = squeeze(nanmean(spike_counts(wake,:),1));
    all_spikes{p,3} = squeeze(nanmean(spike_counts(sleep,:),1));

    all_pearson{p,1} = squeeze(nanmean(pc,1));
    all_pearson{p,2} = squeeze(nanmean(pc(wake,:,:),1));
    all_pearson{p,3} = squeeze(nanmean(pc(sleep,:,:),1));

    all_rl{p,1} = squeeze(nanmean(rl,1));
    all_rl{p,2} = squeeze(nanmean(rl(wake,:),1));
    all_rl{p,3} = squeeze(nanmean(rl(sleep,:),1));

    all_bp{p,1} = squeeze(nanmean(bp,1));
    all_bp{p,2} = squeeze(nanmean(bp(wake,:,:),1));
    all_bp{p,3} = squeeze(nanmean(bp(sleep,:,:),1));

end

%% Save
out.all_names = all_names;
out.all_labels = all_labels;
out.all_coh = all_coh;
out.all_spikes = all_spikes;
out.all_rl = all_rl;
out.all_pearson = all_pearson;
out.all_bp = all_bp;

save([out_folder,'mt_out.mat'],'out')

end