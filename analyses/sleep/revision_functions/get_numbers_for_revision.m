function get_numbers_for_revision

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/seizure_times/'];
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

%% Prep array of numbers
n_one_minute_segments = nan(npts,1);
n_spikes = nan(npts,1);
n_spikes_per_elec = nan(npts,1);
n_seizures = nan(npts,1);
total_dur = nan(npts,1);

for p = 1:npts
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;

    spikes = summ.spikes;
    sz_times = summ.sz_times;
    times = summ.times;

    n_one_minute_segments(p) = size(spikes,2);

    % do the times more or less align
    dur = times(end)-times(1);
    alt_dur = size(spikes,2)*10*60;
    rel_diff = abs(dur-alt_dur)/alt_dur;
    assert(rel_diff < 0.01);
    total_dur(p) = size(spikes,2)/6/24; % in days

    n_spikes(p) = nansum(spikes,'all');
    n_spikes_per_elec(p) = nanmean(nansum(spikes,2));


    n_seizures(p) = size(sz_times,1);


end

error('check it out')

end