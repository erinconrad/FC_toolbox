function sleep_histogram_analysis

%% Parameters
n_periods = 6*3; % want a minimum of 3 hours (18 10-minute blocks) of mostly wake or mostly sleep
min_same = 0.8; % want 80% consistency within that time period
later_search = 12;

time_to_take_spikes = 6*12;
nbins = time_to_take_spikes*2;

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

%% INitialize stuff
all_pts_spikes_bins = nan(npts,nbins);
all_pts_sleep_bins = nan(npts,nbins);
missing_loc = zeros(npts,1);  
rate_anatomy = nan(length(main_locs),npts,nbins);
rate_soz = nan(2,npts,nbins);

for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get main things
    labels = summ.labels;
    ad = summ.ad;
    spikes = summ.spikes;
    loc = summ.ana_loc;
    
    %% Get features for soz vs not
    soz = summ.soz.chs;
    chnums = 1:length(labels);
    is_soz = ismember(chnums,soz);
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    spikes = spikes(~ekg,:);
    loc = loc(~ekg,:);
    is_soz = is_soz(~ekg);
    
    avg_spikes = nanmean(spikes,1);
    
    %% Determine "wake" and "sleep" times
    % normalized ad
    ad_norm = (ad - nanmedian(ad))./iqr(ad);
    wake = ad_norm > disc;
    sleep = ad_norm <= disc;
    
    %% Find transition points and bins
    [transitions,bins] = designate_histogram(sleep,n_periods,min_same,later_search,time_to_take_spikes);
    spikes_in_bins = nan(size(bins));
    sleep_in_bins = nan(size(bins));
    for t = 1:length(transitions)
        spikes_in_bins(t,:) = avg_spikes(bins(t,:));
        sleep_in_bins(t,:) = sleep(bins(t,:));
    end
    
    % Average across transitions
    spikes_in_bins = nanmean(spikes_in_bins,1);
    sleep_in_bins = nanmean(sleep_in_bins,1);
    
    % Add to list of patients
    all_pts_spikes_bins(p,:) = spikes_in_bins;
    all_pts_sleep_bins(p,:) = sleep_in_bins;
    
    %% Bin all spikes
    all_s_bin = nan(size(spikes,1),size(bins,1),size(bins,2));
    for t = 1:length(transitions)
        all_s_bin(:,t,:) = spikes(:,bins(t,:));
    end
    
    %% SOZ vs not
    rate_soz(1,p,:) = nanmean(all_s_bin(is_soz,:,:),[1 2]);
    rate_soz(2,p,:) = nanmean(all_s_bin(~is_soz,:,:),[1 2]);
    
    
    %% Anatomy
    % Skip subsequent loc analyses if missing
    if sum(cellfun(@(x) isempty(x),loc)) == length(loc) 
        missing_loc(p) = 1;
        continue
    end
    
    for sg = 1:length(main_locs)
        ic = ismember(loc,main_locs(sg));
        
       
        % average over anatomical region and transitions
        sp_bins = squeeze(nanmean(all_s_bin(ic,:,:),[1 2]));
        
        % spike rate for that region
        rate_anatomy(sg,p,:) = sp_bins;
    end
    
    
end

figure
set(gcf,'position',[10 10 1400 450])
tiledlayout(1,3,'padding','tight','tilespacing','tight')

%% Main
nexttile
spike_bins = nanmean(all_pts_spikes_bins,1);
sleep_bins = nanmean(all_pts_sleep_bins,1);
times = linspace(-12,12,length(spike_bins));
spp = plot(times,spike_bins,'k-');
hold on
slp = plot(times,sleep_bins,'k--');
title('Spike rate surrounding sleep onset')
xlabel('Hours')
ylabel('Spikes/elecs/min')
xlim([-12 12])
yl = ylim;
plot([0 0],yl,'r--')
legend([spp;slp],{'Spikes','Proportion in sleep'})
ylim(yl);

%% Anatomy
nexttile
anap = nan(size(rate_anatomy,1),1);
for sg = 1:size(rate_anatomy,1)
    anap(sg) = plot(times,squeeze(nanmean(rate_anatomy(sg,:,:),2)));
    hold on
end
title('Spike rate surrounding sleep onset')
xlabel('Hours')
ylabel('Spikes/elecs/min')
xlim([-12 12])
yl = ylim;
plot([0 0],yl,'r--')
legend(anap,main_locs)
ylim(yl);

%% SOZ
nexttile
anap = nan(2,1);
for sg = 1:2
    anap(sg) = plot(times,squeeze(nanmean(rate_soz(sg,:,:),2)));
    hold on
end
title('Spike rate surrounding sleep onset')
xlabel('Hours')
ylabel('Spikes/elecs/min')
xlim([-12 12])
yl = ylim;
plot([0 0],yl,'r--')
legend(anap,main_soz)
ylim(yl);

%% Save
print([out_folder,'histogram'],'-dpng')



end