function out = basic_seeg_analyses

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
edf_out_dir = [results_folder,'edf_out/'];
int_folder = [results_folder,'analysis/backup_intermediate_Feb26_good_spikes/'];
seeg_folder = [results_folder,'analysis/sleep/sleep_seeg/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% possible sleep stages
all_ss = {'R','W','N1','N2','N3'};
nstages = length(all_ss);

%% Initialize variables
rate_ss = nan(npts,nstages);
nseq_ss = nan(npts,nstages);
spread_ss = nan(npts,nstages);
fc_ss = nan(npts,nstages);
all_elec_rates_ss = cell(npts,1);
all_ss = cell(npts,5); %{'R','W','N1','N2','N3'}

%% Loop over patients
for p = 1:npts

    pt_name = strrep(listing(p).name,'.mat','');

    if ~exist([edf_out_dir,pt_name,'/sleep_stage.mat'],'file')
        fprintf('\nskipping %s\n',pt_name);
        continue
    end

    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load spike file
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;

    % Get main things
    spikes = summ.spikes;
    labels = summ.labels;
    ns = summ.ns;
    seq_info = summ.seq_info;
    times = summ.times;
    file_index = summ.file_index;

    %% Load SEEG file
    sout = load([edf_out_dir,pt_name,'/sleep_stage.mat']);
    sout = sout.sout;

    % Get times of sleep transitions
    st_dates= sout.Summary(2:end,2);
    st_times = sout.Summary(2:end,3);
    stage = sout.Summary(2:end,4);
    seeg_secs = convert_transitions_to_ieeg_secs(st_dates,st_times);

    %% Get times with SEEG info and corresponding transitions
    [is_seeg_time,seeg_stage] = match_seeg_times(times,file_index,seeg_secs,stage);

    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    spikes = spikes(~ekg,:);
    labels = labels(~ekg);
    ns = ns(~ekg,:);

    %% Spike rate in different sleep stages

    % take average of spikes across electrodes
    mean_spikes = nanmean(spikes,1); %spikes/elecs/min
    
    % spike rate across sleep stages
    for is = 1:nstages
        rate_ss(p,is) = nanmean(mean_spikes(strcmp(seeg_stage,all_ss{is})));
    end

    %% Spike sequence info and ns info
    mean_ns = squeeze(nanmean(ns,1)); % node strength averaged across electrodes
    for is = 1:nstages
        nseq_ss(p,is) = nanmean(seq_info(1,strcmp(seeg_stage,all_ss{is})));
        spread_ss(p,is) = nanmean(seq_info(2,strcmp(seeg_stage,all_ss{is})));
        fc_ss(p,is) = nanmean(mean_ns(strcmp(seeg_stage,all_ss{is})));
    end

    %% All electrode spike rates in different sleep stages (for model)
    for is = 1:nstages
        all_elec_rates_ss{p}(:,is) = nanmean(spikes(:,strcmp(seeg_stage,all_ss{is})),2);
        all_ss{p,is} = strcmp(seeg_stage,all_ss{is});
    end

end

out.all_elec_rates_ss = all_elec_rates_ss;
out.fc_ss = fc_ss;
out.spread_ss = spread_ss;
out.nseq_ss = nseq_ss;
out.rate_ss = rate_ss;

end