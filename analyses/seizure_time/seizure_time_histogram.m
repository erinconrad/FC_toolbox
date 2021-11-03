function seizure_time_histogram(rm_cluster)

%% Parameters
surround_hours = 12;


surround_secs = surround_hours*3600;
surround = surround_hours*6;
nbins = surround*2;

if rm_cluster == 1
    rm_cluster_text = '_rm_sz_clusters';
else
    rm_cluster_text = '_keep_sz_clusters';
end


%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/seizure_times/'];
int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% initialize
all_pts_spikes_bins = nan(npts,nbins);

for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get main things
    labels = summ.labels;
    spikes = summ.spikes;
    sz_times = summ.sz_times;
    times = summ.times;
    times = times;
    sz_times = sz_times;
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    
    spikes = spikes(~ekg,:);
    avg_spikes = nanmean(spikes,1);
    
    
    %% test plots
    if 0
        figure
        h = turn_nans_gray(spikes);
        hold on
        set(h,'XData',times);
        for s = 1:size(sz_times,1)
            plot([sz_times(s,1) sz_times(s,1)],ylim,'r--')
            plot([sz_times(s,2) sz_times(s,2)],ylim,'r--')
        end
        xlim([times(1) times(end)])
        
        figure
        plot(times,avg_spikes)
        hold on
        for s = 1:size(sz_times,1)
            plot([sz_times(s,1) sz_times(s,1)],ylim,'r--')
            plot([sz_times(s,2) sz_times(s,2)],ylim,'r--')
        end
        xlim([times(1) times(end)])
    end
    
    %% bin times around seizures
    bins = bin_seizure_times(sz_times,times,surround,rm_cluster,surround_secs);
    
    %% Get spikes in those times
    sp_bins = nan(size(bins));
    for s = 1:size(bins,1)
        
        curr_bins = bins(s,:);
        
        % get only those bins that are within range
        bins_in_range = curr_bins > 0 & curr_bins <= size(avg_spikes,2);
        too_early = curr_bins <= 0;
        too_late = curr_bins > size(avg_spikes,2);
        
        % get the spikes for those in range
        out_bins = [nan(1,sum(too_early)),avg_spikes(bins_in_range),nan(1,sum(too_late))];
        
        sp_bins(s,:) = out_bins;
    end
    
    %% Average over seizures
    sz_avg_sp_bins = nanmean(sp_bins,1);
    
    %% Add to patient matrix
    all_pts_spikes_bins(p,:) = sz_avg_sp_bins;
end

figure
sp_bins = nanmean(all_pts_spikes_bins,1);
times = linspace(-surround_hours,surround_hours,length(sp_bins));
plot(times,sp_bins,'k');
hold on
title('Spike rate surrounding seizure onset')
xlabel('Hours')
ylabel('Spikes/elecs/min')
xlim([-surround_hours surround_hours])
yl = ylim;
plot([0 0],yl,'r--')
ylim(yl)
print([out_folder,'histogram',rm_cluster_text],'-dpng')
close(gcf)

end