function out = seizure_time_histogram(disc)

%{
Build a model:
spike rate ~ time relative to seizure + sleep vs wake + (1|patient)
   - AR???
   - logistic??
%}

%% Parameters
surround_hours = 6;
surround_secs = surround_hours*3600; % convert to seconds
surround = surround_hours*6; % convert to #bins (6 hours = 6*6 10 minute bins)
nbins = surround*2;
rm_cluster = 0;

main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
main_lats = {'Left','Right'};
main_soz = {'SOZ','Not SOZ'};
main{1} = main_locs;
main{2} = main_lats;
main{3} = main_soz;


if rm_cluster == 1
    rm_cluster_text = '_rm_clust';
else
    rm_cluster_text = '_keep_clust';
end

%{
if do_avg == 1
    do_avg_text = '_avg';
else
    do_avg_text = '_no_avg';
end
%}


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
all_pts_sleep_bins = nan(npts,nbins);

all_pts_sp_vec = [];
all_pts_sz_id_vec = [];
all_pts_bin_id_vec = [];
all_pts_id_vec = [];
all_pts_sleep_vec = [];
missing_loc = zeros(npts,1);
post_ictal_rates = cell(npts,1);
pre_to_post_change = nan(npts,1);

spikes_strat = cell(3,1);
for i = 1:length(spikes_strat)
    spikes_strat{i} = nan(length(main{i}),npts,nbins);
end

names = cell(npts,1);


% start running count of which seiuzre
sz_count = 0;
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
    loc = summ.ana_loc;
    lat = summ.ana_lat;
    name = summ.name;
    
    names{p} = name;
    
    % Fix lat thing
    for i = 1:length(lat)
        if isempty(lat{i}), lat{i} = 'unspecified'; end
    end
    
    %% Get features for soz vs not
    soz = summ.soz.chs;
    chnums = 1:length(labels);
    is_soz = ismember(chnums,soz);
    
    %% Find and remove non-intracranial
    ekg = find_non_intracranial(labels);
    
    spikes = spikes(~ekg,:);
    avg_spikes = nanmean(spikes,1);
    spike_counts = nansum(spikes,1);
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    loc = loc(~ekg,:);
    lat = lat(~ekg);
    is_soz = is_soz(~ekg);
    soz_text = cell(sum(~ekg),1);
    soz_text(is_soz) = {'SOZ'};
    soz_text(~is_soz) = {'Not SOZ'};
    
    %% Determine "wake" and "sleep" times
    % normalized ad
    ad_norm = (ad - nanmedian(ad))./iqr(ad);
    sleep = ad_norm <= disc;
    sleep(isnan(sleep)) = 0;
    
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
    sp_count_bins = nan(size(bins));
    sp_count_bins_all_elecs = nan(size(bins,1),size(spikes,1));
    
    for s = 1:size(bins,1)
        
        curr_bins = bins(s,:);
        
        % get only those bins that are within range
        bins_in_range = curr_bins > 0 & curr_bins <= size(avg_spikes,2);
        too_early = curr_bins <= 0;
        too_late = curr_bins > size(avg_spikes,2);
        
        % get the spikes for those in range
        out_bins = [nan(1,sum(too_early)),avg_spikes(curr_bins(bins_in_range)),nan(1,sum(too_late))];
        sp_bins(s,:) = out_bins;
        
        sp_count_bins(s,:) = [nan(1,sum(too_early)),spike_counts(curr_bins(bins_in_range)),nan(1,sum(too_late))];
        sp_count_bins_all_elecs(s,:) = nanmean(spikes(:,curr_bins(bins_in_range)),2);
        
        % sleep bins
        sleep_in_bins(s,:) = [nan(1,sum(too_early)),...
            sleep(curr_bins(bins_in_range)),nan(1,sum(too_late))];

    end
    
    post_ictal_rates{p} = (nanmean(sp_count_bins_all_elecs,1))';
    
    %% Calculate the pre-to-post change
    pre_bins = 1:nbins/2;
    post_bins = nbins/2+1:nbins;
    pre_to_post_change(p) = nanmean(sp_count_bins(:,post_bins) - sp_count_bins(:,pre_bins),'all');
    
    %% Get spike rate in each bin for each group for locs and lats
    % Loop over loc vs lat vs soz
    for g = 1:3
        if g == 1
            group = loc;
            % Skip subsequent loc analyses if missing
            if sum(cellfun(@(x) isempty(x),loc)) == length(loc) 
                missing_loc(p) = 1;
                continue
            end
        elseif g == 2
            group = lat;
        elseif g == 3
            group = soz_text;
        end
        
        % Get the rates corresponding to the subgroups
        % (can probably do this without a for loop)
       
        for sg = 1:length(main{g})
            ic = ismember(group,main{g}(sg));
            
            sp_bins_loc = nan(size(bins));
            
            for s = 1:size(bins,1)
                
                curr_bins = bins(s,:);
        
                % get only those bins that are within range
                bins_in_range = curr_bins > 0 & curr_bins <= size(avg_spikes,2);
                too_early = curr_bins <= 0;
                too_late = curr_bins > size(avg_spikes,2);
                
                % get the spikes for those in range for correct locs
                sp_bins_loc(s,:) = [nan(1,sum(too_early)),...
                    nanmean(spikes(ic,curr_bins(bins_in_range)),1),nan(1,sum(too_late))];
        
            end
            
            % Average across seizures
            sp_bins_loc = nanmean(sp_bins_loc,1);
            
            % add to cell array
            spikes_strat{g}(sg,p,:) =  sp_bins_loc;
            
        end
        

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
    
    %% Prep for model
    nszs = size(sp_bins,1);
    sz_ids = (1:nszs)' + sz_count;
    sz_count = sz_count + nszs; % increase total seizure count (for next patient)
    sz_id_bins = repmat(sz_ids,1,nbins);
    bin_id_bins = repmat(1:nbins,nszs,1);
    
    % reshape things to one dimensional vector
    sp_vec = reshape(sp_count_bins,nszs*nbins,1); % spike rates
    sz_id_vec = reshape(sz_id_bins,nszs*nbins,1); % identify which seizure
    bin_id_vec = reshape(bin_id_bins,nszs*nbins,1); % which bin (will treat as categorical)
    sleep_vec = reshape(sleep_in_bins,nszs*nbins,1); % sleep or wake (binary)
    
    % patient identifier
    pt_id_vec = p*ones(length(sp_vec),1);
    
    % Add these to overall pt
    all_pts_sp_vec = [all_pts_sp_vec;sp_vec];
    all_pts_sz_id_vec = [all_pts_sz_id_vec;sz_id_vec];
    all_pts_bin_id_vec = [all_pts_bin_id_vec;bin_id_vec];
    all_pts_id_vec = [all_pts_id_vec;pt_id_vec];
    all_pts_sleep_vec = [all_pts_sleep_vec;sleep_vec];
    
    %% Average over seizures
    sz_avg_sp_bins = nanmean(sp_bins,1);
    sleep_in_bins = nanmean(sleep_in_bins,1);
    
    
    %% Add to patient matrix
    all_pts_spikes_bins(p,:) = sz_avg_sp_bins;
    all_pts_sleep_bins(p,:) = sleep_in_bins;
end

%% Remove empty pts for loc analyses
missing_loc = logical(missing_loc);
spikes_strat{1}(:,missing_loc,:) = [];


%% Do model
%{

if do_avg
    % Make bin identifier
    bin_id = repmat(1:nbins,npts,1);
    
    % Make pt identifier
    pt_id = repmat((1:npts)',1,nbins);
    
    % Reshape into one dimensional vector
    npts = size(all_pts_spikes_bins,1);
    spikes_vec = reshape(all_pts_spikes_bins,npts*nbins,1);
    sleep_vec = reshape(all_pts_sleep_bins,npts*nbins,1);
    bin_vec = reshape(bin_id,npts*nbins,1);
    pt_vec = reshape(pt_id,npts*nbins,1);
    
    T = table(spikes_vec,bin_vec,sleep_vec,pt_vec,...
        'VariableNames',...
        {'SpikeRate','Bin','Sleep','Patient'});
    
    % Establish categorical variables (note that sleep is now
    % pseudo-linear, not categorical).
    T.Bin = categorical(T.Bin);
    T.Patient = categorical(T.Patient);
    
    lme = fitlme(T,'SpikeRate ~ Bin + Sleep + (1|Patient)');
    
else
    % Treat every seizure separately
    T = table(all_pts_sp_vec,all_pts_bin_id_vec,(all_pts_sleep_vec),...
        all_pts_sz_id_vec,all_pts_id_vec,'VariableNames',...
        {'SpikeRate','Bin','Sleep','Seizure','Patient'});

    % establish categorical variables
    T.Bin = categorical(T.Bin);
    T.Sleep = categorical(T.Sleep);
    T.Seizure = categorical(T.Seizure);
    T.Patient = categorical(T.Patient);


    lme = fitlme(T,'SpikeRate~ Bin + Sleep + (1|Seizure) + (1|Patient)');

    lme_no_sleep = fitlme(T,'SpikeRate~ Bin + (1|Seizure) + (1|Patient)');
end

% Find the ids of the significant bins (relative to the first bin)
sig_bins = lme.Coefficients.pValue < 0.05/nbins; % Bonferroni correction
sig_bins(1) = 0; % the first one is the intercept, ignore
sig_bins = sig_bins(1:end-1); % the last one is the coefficient for sleep, remove
sp_bins = nanmean(all_pts_spikes_bins,1);
%}
times = linspace(-surround_hours,surround_hours,length(sp_bins));
%out.lme = lme;
%out.sig_bins = sig_bins;
out.all_pts_spikes_bins = all_pts_spikes_bins;
out.all_pts_sleep_bins = all_pts_sleep_bins;
out.times = times;
out.surround_hours = surround_hours;
%out.xlim = [-surround_hours surround_hours];
%out.T = T;
out.spikes_strat = spikes_strat;
out.main = main;
out.names = names;
out.post_ictal_rates = post_ictal_rates;
out.pre_to_post_change = pre_to_post_change;

%{
figure
tiledlayout(2,1,'tilespacing','tight','padding','tight')

nexttile
sp_bins = nanmean(all_pts_spikes_bins,1);
sleep_bins = nanmean(all_pts_sleep_bins,1);
times = linspace(-surround_hours,surround_hours,length(sp_bins));
psp = plot(times,sp_bins);
hold on
plot(times(sig_bins),sp_bins(sig_bins)+0.05,'r*');
title('Spike rate surrounding seizure onset')
ylabel('Spikes/min')
xlabel('Hours')
xlim([-surround_hours surround_hours])
plot([0 0],ylim,'r--')

nexttile
psl = plot(times,sleep_bins);
hold on
title('Proportion detected to be asleep surrounding seizure onset')
xlabel('Hours')
ylabel('Proportion detected to be asleep')
xlim([-surround_hours surround_hours])
plot([0 0],ylim,'r--')

print([out_folder,'histogram',rm_cluster_text,do_avg_text],'-dpng')
close(gcf)
%}

end