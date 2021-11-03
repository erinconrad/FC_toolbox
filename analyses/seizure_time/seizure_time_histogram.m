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

%% AD validation
[roc,auc,disc] = ad_validation;

%% initialize
all_pts_spikes_bins = nan(npts,nbins);
all_pts_sleep_bins = nan(npts,nbins);

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
    ad = summ.ad;
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    
    spikes = spikes(~ekg,:);
    avg_spikes = nanmean(spikes,1);
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    
    %% Determine "wake" and "sleep" times
    % normalized ad
    ad_norm = (ad - nanmedian(ad))./iqr(ad);
    wake = ad_norm > disc;
    sleep = ad_norm <= disc;
    
    
    %% bin times around seizures
    bins = bin_seizure_times(sz_times,times,surround,rm_cluster,surround_secs);
    
    %% test plots
    if 0
        
        figure
        plot(times,avg_spikes)
        hold on
        for s = 1:size(sz_times,1)
            plot([sz_times(s,1) sz_times(s,1)],ylim,'b-','linewidth',2)
        end
        for s = 1:size(bins,1)
            plot([times(bins(s,1)) times(bins(s,1))],ylim,'g--','linewidth',1)
            plot([times(bins(s,end)) times(bins(s,end))],ylim,'r--','linewidth',1)
        end
        xlim([times(1) times(end)])
    end
    
    %% Get spikes and sleep in those times
    
    sp_bins = nan(size(bins));
    sleep_in_bins = nan(size(bins));

    for s = 1:size(bins,1)
        
        curr_bins = bins(s,:);
        
        % get only those bins that are within range
        bins_in_range = curr_bins > 0 & curr_bins <= size(avg_spikes,2);
        too_early = curr_bins <= 0;
        too_late = curr_bins > size(avg_spikes,2);
        
        % get the spikes for those in range
        out_bins = [nan(1,sum(too_early)),avg_spikes(curr_bins(bins_in_range)),nan(1,sum(too_late))];
        sp_bins(s,:) = out_bins;
        
        
        % sleep bins
        sleep_in_bins(s,:) = [nan(1,sum(too_early)),...
            sleep_in_bins(curr_bins(bins_in_range)),nan(1,sum(too_late))];

    end
    
    %% test plots
    if 0
        
        figure
        plot(times,avg_spikes)
        hold on
        for s = 1:size(sz_times,1)
            plot([sz_times(s,1) sz_times(s,1)],ylim,'b-','linewidth',2)
        end
        for s = 1:size(bins,1)
            plot(times(bins(s,:)),sp_bins(s,:),'ro');
            %plot(times(bins(s,:)),avg_spikes(bins(s,:)),'ro');
        end
        xlim([times(1) times(end)])
    end
    
    %% Average over seizures
    sz_avg_sp_bins = nanmean(sp_bins,1);
    sleep_in_bins = nanmean(sleep_in_bins,1);
    
    
    %% Add to patient matrix
    all_pts_spikes_bins(p,:) = sz_avg_sp_bins;
    all_pts_sleep_bins(p,:) = sleep_in_bins;
end

figure
sp_bins = nanmean(all_pts_spikes_bins,1);
sleep_bins = nanmean(all_pts_sleep_bins,1);
times = linspace(-surround_hours,surround_hours,length(sp_bins));
psp = plot(times,sp_bins);
hold on
psl = plot(times,sleep_bins);
title('Spike rate surrounding seizure onset')
xlabel('Hours')
xlim([-surround_hours surround_hours])
yl = ylim;
plot([0 0],yl,'r--')
legend([psp,psl],{'Spikes/elec/min','Proportion asleep'});
ylim(yl)
print([out_folder,'histogram',rm_cluster_text],'-dpng')
close(gcf)

end