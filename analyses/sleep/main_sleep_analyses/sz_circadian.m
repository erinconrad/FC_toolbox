function out = sz_circadian(disc,exc)

%{
This analysis obtains seizure counts as a function of time of day
%}

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
int_folder = [results_folder,'analysis/intermediate/'];
out_folder = [results_folder,'analysis/sleep/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% Get the longest run (will pad the others with zeros)
longest_run = 0;
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    run_length = length(summ.times);
    if run_length > longest_run
        longest_run = run_length;
    end
end

%% Initialize PSD
all_psd = nan(npts,ceil(longest_run/2));
fs = 0.0017;%1/summ(1).block_dur;
all_freqs = nan(npts,ceil(longest_run/2));
all_pre_wake = cell(npts,1);
sz_rate_sw = nan(npts,2);


%% for mod midnight, get file size
[~,n_tod_bins,tod_edges] = bin_mod_midnight_times(zeros(5000,1),[]);
all_tod_rate = nan(npts,n_tod_bins); %w, s
    
%% Loop over patients and get psd per pt
for p = 1:npts
    
    %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    times = summ.times;
    run_length = length(times);
    
    if isempty(summ.sz_times)
        continue;
        
    end
    
    sz_times = summ.sz_times(:,1);
    block = summ.block_dur;
    %time_bins = block:block:ceil(times(end)/block)*block;
    ad = summ.ad;
    labels = summ.labels;
    mod_midnight = summ.mod_midnight;
    
    
    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    
    %% Determine wake and sleep
    [sleep,wake] = find_sleep_wake(ad,exc,disc);
    
    %% Bin sz times to get counts per bin (in the original time bins)
    time_bins = times;
    for t = 1:length(times)
        time_bins(t) = max(times(t),ceil(times(t)/block)*block);
    end
    [bin_counts,preceding_bin_counts,closest_bins] = bin_stuff(sz_times,time_bins);
    
    %% Get sz rate by time of day
     % Bin the mod midnights
    [mod_midnight,nbins,edges] = bin_mod_midnight_times(mod_midnight,tod_edges);
    tod_rate = nan(n_tod_bins,1);
    for t = 1:n_tod_bins
        curr_bins = mod_midnight == t; % which runs match that time of day
        tod_rate(t,:) = nansum(bin_counts(curr_bins),'all'); % how many seizures are in that time of day bin

            
    end
    all_tod_rate(p,:) = tod_rate;
    
    %% Get sz frequency in wake and sleep
    sz_rate_sw(p,:) = [nanmean(preceding_bin_counts(wake)) nanmean(preceding_bin_counts(sleep))];
    
    %% Preceding bin wake vs sleep
    closest_bins = closest_bins(~isnan(closest_bins));
    pre_bins = closest_bins - 1;
    pre_wake = nan(length(closest_bins),1);
    pre_wake(wake(pre_bins)==1)=1;
    pre_wake(sleep(pre_bins)==1)=0;
    all_pre_wake{p} = pre_wake;
    
    %% pad spikes for PSD
    bin_counts = [bin_counts;zeros(longest_run-run_length,1)];
    
    %% get psd
    [P,freqs] = power_by_freq(bin_counts,fs);
    
    all_psd(p,:) = P;
    all_freqs(p,:) = freqs;
    
end


%% Stuff
periods = 1./freqs/3600;
low_period = periods <= 100;
periods = periods(low_period);

% Restrict to less than 100
% Get stats
all_psd = all_psd(:,low_period);
all_psd = all_psd./sum(all_psd,2);
median_psd = median(all_psd,1);
iqr_psd = [prctile(all_psd,25,1);prctile(all_psd,75,1)];

out.all_psd = all_psd;
out.low_period = low_period;
out.median_psd = median_psd;
out.iqr_psd = iqr_psd;
out.periods = periods;
out.pre_wake = all_pre_wake;
out.sz_rate_sw = sz_rate_sw;
out.all_tod_rate = all_tod_rate;

%mp = shaded_error_bars(periods,median_psd,iqr_psd,[0 0 0]);

end