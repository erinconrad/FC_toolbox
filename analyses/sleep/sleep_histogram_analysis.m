function out = sleep_histogram_analysis(rm_cluster,disc)

%{
I need to do a linear rather than a poisson model because I don't have
counts perse but sums of counts across channels. If I really wanted to I
could go back at the level of grouping spikes in bins and define "meta
spike counts" to be a "meta spike" anytime there is a spike in any channel
for a defined time period. But the linear model seems ok.
%}

%% Parameters
n_periods = 6*3; % want a minimum of 3 hours (18 10-minute blocks) of mostly wake or mostly sleep
min_same = 0.7; % want 70% consistency within that time period
later_search = 12;

time_to_take_spikes = 6*12;
nbins = time_to_take_spikes*2;

main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
main_lats = {'Left','Right'};
main_soz = {'SOZ','Not SOZ'};
main{1} = main_locs;
main{2} = main_lats;


if rm_cluster == 1
    rm_cluster_text = '_rm_clusters';
else
    rm_cluster_text = '_keep_clusters';
end

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

ad_plot_folder = [out_folder,'sleep_disc/'];
if ~exist(ad_plot_folder,'dir'), mkdir(ad_plot_folder); end

%% Load summary file
%{
summ = load([summ_folder,'summ.mat']);
summ = summ.summ;
%}

%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% INitialize stuff
all_pts_spikes_bins = nan(npts,nbins);
all_pts_sleep_bins = nan(npts,nbins);
missing_loc = zeros(npts,1);  
rate_anatomy = nan(length(main_locs),npts,nbins);
rate_soz = nan(2,npts,nbins);
pt_sz_bins = cell(npts,1);
all_sz_bins = [];
all_transitions = cell(npts,1);

all_pts_sp_vec = [];
all_pts_bin_id_vec = [];
all_pts_id_vec = [];
all_pts_trans_id_vec = [];

names = cell(npts,1);

% start running count of which sleep transition
trans_count = 0;
for p = 1:npts
    
    
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get main things
    labels = summ.labels;
    ad = summ.ad;
    spikes = summ.spikes;
    loc = summ.ana_loc;
    name = summ.name;
    sz_times = summ.sz_times;
    times = summ.times;
    
    names{p} = name;
    
    fprintf('\nDoing %s, patient %d of %d\n',name,p,npts);
    
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
    sp_counts = nansum(spikes,1);
    
    %% Determine "wake" and "sleep" times
    % normalized ad
    ad_norm = (ad - nanmedian(ad))./iqr(ad);
    wake = ad_norm > disc;
    sleep = ad_norm <= disc;
    
    %% Find transition points and bins
    [transitions,bins] = designate_histogram(sleep,wake,n_periods,min_same,...
        later_search,time_to_take_spikes,rm_cluster,ad_norm,disc,name,ad_plot_folder);
    spikes_in_bins = nan(size(bins));
    sleep_in_bins = nan(size(bins));
    sp_counts_in_bins = nan(size(bins));
    for t = 1:length(transitions)
        spikes_in_bins(t,:) = avg_spikes(bins(t,:));
        sleep_in_bins(t,:) = sleep(bins(t,:));
        sp_counts_in_bins(t,:) = sp_counts(bins(t,:));
    end
    all_transitions{p} = transitions;
    
    % Average across transitions
    spikes_in_bins = nanmean(spikes_in_bins,1);
    sleep_in_bins = nanmean(sleep_in_bins,1);
    
    % Add to list of patients
    all_pts_spikes_bins(p,:) = spikes_in_bins;
    all_pts_sleep_bins(p,:) = sleep_in_bins;
    
    %% bin seizure times into this framework
    if ~isempty(sz_times)
        sz_start = sz_times(:,1);
        sz_bins = find_segment_closest_time(bins,times,sz_start);

        pt_sz_bins{p} = sz_bins;
        all_sz_bins = [all_sz_bins;sz_bins];
    end
    
    %% Also prep for a model
    nbins = size(sp_counts_in_bins,2);
    ntrans = size(sp_counts_in_bins,1);
    trans_ids = (1:ntrans)' + trans_count;
    trans_count = trans_count + ntrans; % increase total sleep transition count (for next patient)
    trans_id_bins = repmat(trans_ids,1,nbins);
    bin_id_bins = repmat(1:nbins,ntrans,1);
    
    % reshape things to one dimensional vector
    sp_vec = reshape(sp_counts_in_bins,ntrans*nbins,1); % spike rates
    trans_id_vec = reshape(trans_id_bins,ntrans*nbins,1); % identify which sleep transition
    bin_id_vec = reshape(bin_id_bins,ntrans*nbins,1); % which bin (will treat as categorical)
    pt_id_vec = p*ones(length(sp_vec),1);
    
    % Add to overall pt
    all_pts_sp_vec = [all_pts_sp_vec;sp_vec];
    all_pts_bin_id_vec = [all_pts_bin_id_vec;bin_id_vec];
    all_pts_id_vec = [all_pts_id_vec;pt_id_vec];
    all_pts_trans_id_vec = [all_pts_trans_id_vec;trans_id_vec];
    
    
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

%% Model
T = table(all_pts_sp_vec,all_pts_bin_id_vec,all_pts_trans_id_vec,all_pts_id_vec,...
    'VariableNames',{'SpikeRate','Bin','SleepTransition','Patient'});

% establish categorical variables
T.Bin = categorical(T.Bin);
T.SleepTransition = categorical(T.SleepTransition);
T.Patient = categorical(T.Patient);

% Do a linear mixed effects model, treating each bin as a separate
% categorical predictor
lme = fitlme(T,'SpikeRate~ Bin + (1|SleepTransition) + (1|Patient)');

% Find the ids of the significant bins (relative to the first bin)
sig_bins = lme.Coefficients.pValue < 0.05/nbins; % Bonferroni correction
sig_bins(1) = 0; % the first one is the intercept, ignore

spike_bins = nanmean(all_pts_spikes_bins,1);
times = linspace(-12,12,length(spike_bins));
out.lme = lme;
out.sig_bins = sig_bins;
out.all_pts_spikes_bins = all_pts_spikes_bins;
out.all_pts_sleep_bins = all_pts_sleep_bins;
out.times = times;
out.title = 'Spike rate surrounding sleep onset';
out.xlabel = 'Hours';
out.ylabel = 'Spikes/elecs/min';
out.xlim = [-12 12];
out.T = T;
out.names = names;
out.pt_sz_bins = pt_sz_bins;
out.all_sz_bins = all_sz_bins;
out.all_transitions = all_transitions;

%{
lme_with_pt = fitlme(T,'SpikeRate~ Bin + (1|SleepTransition) + (1|Patient)');
lme_with_pt
%}
%{
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
plot(times(sig_bins),spike_bins(sig_bins)+0.05,'r*');
%hold on
%slp = plot(times,sleep_bins,'k--');
title('Spike rate surrounding sleep onset')
xlabel('Hours')
ylabel('Spikes/elecs/min')
xlim([-12 12])
yl = ylim;
plot([0 0],yl,'r--')
%legend([spp;slp],{'Spikes','Proportion in sleep'})
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
print([out_folder,'histogram',rm_cluster_text],'-dpng')
close(gcf)
%}

end