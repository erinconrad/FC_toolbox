function basic_seeg_analyses

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
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

%% Loop over patients
for p = 1:npts


    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load spike file
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;

    % Get main things
    loc = summ.ana_loc;
    lat = summ.ana_lat;
    spikes = summ.spikes;
    ad = summ.ad;
    rl = summ.rl;
    coi_global = summ.coi_global;
    labels = summ.labels;
    ns = summ.ns;
    name = summ.name;
    seq_info = summ.seq_info;
    leader = summ.leader;
    mod_midnight = summ.mod_midnight;

    %% Load SEEG file
    sout = load([edf_out_dir,pt_name,'/sleep_stage.mat']);
    sout = sout.sout;

    % Get times of sleep transitions
    st_dates= sout.Summary(2:end,2);
    st_times = sout.Summary(2:end,3);
    stage = sout.Summary(2:end,4);
    seeg_secs = convert_transitions_to_ieeg_secs(st_dates,st_times);


end

end