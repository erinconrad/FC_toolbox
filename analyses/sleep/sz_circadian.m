function out = sz_circadian(disc,exc)


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
    
    
    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    
    %% Determine wake and sleep
    [sleep,wake] = find_sleep_wake(ad,exc,disc);
    
    %% Bin sz times
    time_bins = times;
    for t = 1:length(times)
        time_bins(t) = max(times(t),ceil(times(t)/block)*block);
    end
    [bin_counts,preceding_bin_counts,closest_bins] = bin_stuff(sz_times,time_bins);
    
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
    bin_counts = [bin_counts,zeros(1,longest_run-run_length)];
    
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

%mp = shaded_error_bars(periods,median_psd,iqr_psd,[0 0 0]);

end