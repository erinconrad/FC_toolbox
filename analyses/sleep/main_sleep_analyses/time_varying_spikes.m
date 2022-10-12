function out = time_varying_spikes(disc)

exc = []; % don't change

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate_epilepsia_revision/'];%[results_folder,'analysis/backup_intermediate_Feb26_good_spikes/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% Initialize stuff
all_spikes = cell(npts,1);
all_ws = cell(npts,2);
all_is_soz = cell(npts,1);
all_elecs_names = cell(npts,1);
all_names = cell(npts,1);

%% Loop over patients
for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get main things
    spikes = summ.spikes;
    ad = summ.ad;
    labels = summ.labels;
    name = summ.name;
    
    %% Get features for soz vs not
    soz = summ.soz.chs;
    chnums = 1:length(labels);
    is_soz = ismember(chnums,soz);
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    spikes = spikes(~ekg,:);
    is_soz = is_soz(~ekg);
    
    %% Determine "wake" and "sleep" times
    [sleep,wake] = find_sleep_wake(ad,exc,disc);

    %% Fill up
    all_names{p} = name;
    all_spikes{p} = spikes;
    all_is_soz{p} = is_soz;
    all_ws{p,1} = wake;
    all_ws{p,2} = sleep;
    all_elecs_names{p} = labels;
    
    
end

out.all_names = all_names;
out.all_elecs_names = all_elecs_names;
out.all_ws = all_ws;
out.all_spikes = all_spikes;
out.all_is_soz = all_is_soz;

end