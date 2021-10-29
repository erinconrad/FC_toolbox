function sleep_histogram_analysis

%% Parameters
n_periods = 6*3; % want a minimum of 3 hours (18 10-minute blocks) of mostly wake or mostly sleep
min_same = 0.8; % want 80% consistency within that time period
later_search = 12;

time_to_take_spikes = 6*12;
half_cycle = 6*12;

main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
main_lats = {'Left','Right'};
main_soz = {'SOZ','Not SOZ'};
main{1} = main_locs;
main{2} = main_lats;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load summary file
%{
summ = load([summ_folder,'summ.mat']);
summ = summ.summ;
%}

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% AD validation
[roc,auc,disc] = ad_validation;

all_pts_bins = [];

for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    labels = summ.labels;
    ad = summ.ad;
    ekg = find_non_intracranial(labels);
    spikes = summ.spikes;
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    spikes = spikes(~ekg,:);
    
    avg_spikes = nanmean(spikes,1);
    
    %% Determine "wake" and "sleep" times
    % normalized ad
    ad_norm = (ad - nanmedian(ad))./iqr(ad);
    wake = ad_norm > disc;
    sleep = ad_norm <= disc;
    
    %% Find transition points and bins
    [transitions,bins] = designate_histogram(sleep,n_periods,min_same,later_search,time_to_take_spikes);
    spikes_in_bins = nan(size(bins));
    for t = 1:length(transitions)
        spikes_in_bins(t,:) = avg_spikes(bins(t,:));
    end
    
    all_pts_bins = [all_pts_bins;spikes_in_bins];
    
end

if 1
    figure
    bins = nanmean(all_pts_bins,1);
    times = linspace(-12,12,length(bins));
    plot(times,bins)
    hold on
    plot([0 0],ylim,'r--')
    title('Spike rate surrounding sleep onset')
    xlabel('Hours')
    ylabel('Spikes/elecs/min')
    xlim([-12 12])
    print([out_folder,'histogram'],'-dpng')
end

end