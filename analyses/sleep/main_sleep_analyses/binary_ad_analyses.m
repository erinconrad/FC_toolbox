function out = binary_ad_analyses(disc,exc)

%{
This is the function that does most analyses comparing what happens to
spikes in wake and sleep. It takes disc, the normalized alpha delta ratio
used as the cutoff between wake and sleep, and exc, any time points to be
excluded (this should be empty as I am not excluding time points).
%}

%% Parameters
%min_spikes = 0.1;

main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
main_lats = {'Left','Right'};
main_soz = {'SOZ','Not SOZ'};
main{1} = main_locs;
main{2} = main_lats;
main{3} = main_soz;

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


%% Listing of available files
listing = dir([int_folder,'*.mat']);
npts = length(listing);

%% Main analyses
%{
r_ad_ana = cell(3,1);
for i = 1:length(r_ad_ana)
    r_ad_ana{i} = nan(length(main{i}),npts,2);
end

rate_strat_ana = cell(3,1);
for i = 1:length(rate_strat_ana)
    rate_strat_ana{i} = nan(length(main{i}),npts,2);
end


rl_strat_ana = cell(3,1);
for i = 1:length(rl_strat_ana)
    rl_strat_ana{i} = nan(length(main{i}),npts,2);
end

r_rl_ana = cell(3,1);
for i = 1:length(r_rl_ana)
    r_rl_ana{i} = nan(length(main{i}),npts,2);
end

rate_overall_ana = cell(3,1);
for i = 1:length(rate_overall_ana)
    rate_overall_ana{i} = nan(length(main{i}),npts);
end

rl_overall_ana = cell(3,1);
for i = 1:length(rl_overall_ana)
    rl_overall_ana{i} = nan(length(main{i}),npts);
end
%}

%% Initialize other stuff
null_ps = nan(npts,1);
ns_sw = nan(npts,2);
all_rates = nan(npts,2);
overall_rates = nan(npts,1);
all_coi = nan(npts,2);
missing_loc = zeros(npts,1);
n_sleep_wake = nan(npts,2);
names = cell(npts,1);
seq_sw = nan(npts,4);
%soz_rank_sw = nan(npts,2);
%soz_rank_sw_rl = nan(npts,2);
%rl_sw_corr = nan(npts,1);
nspikey = nan(npts,2);
all_is_soz = cell(npts,1);
all_elecs_rates = cell(npts,1);
all_elecs_rl = cell(npts,1);
all_elecs_names = cell(npts,1);

all_elecs_rates_sw = cell(npts,1);
all_elecs_rl_sw = cell(npts,1);

all_elecs_leader = cell(npts,1);
all_elecs_leader_sw = cell(npts,1);

all_elecs_ns_sw = cell(npts,1);

%% for mod midnight, get file size
[~,n_tod_bins,tod_edges] = bin_mod_midnight_times(zeros(5000,1),[]);
all_tod_sw = nan(npts,n_tod_bins,2); %w, s

%% Loop over patients
for p = 1:npts
    
    fprintf('\nDoing patient %d of %d\n',p,npts);
    
     %% Load
    summ = load([int_folder,listing(p).name]);
    summ = summ.summ;
    
    %% Get main things
    loc = summ.ana_loc;
    lat = summ.ana_lat;
    spikes = summ.spikes;
    ad = summ.ad;
    rl = summ.rl;
    coi_global = summ.coi_global;
    labels = summ.labels;
    ns = summ.ns;
    name = summ.name;
    seq_info = summ.seq_info;
    leader = summ.leader;
    mod_midnight = summ.mod_midnight;
    
    names{p} = name;
    
    % Bin the mod midnights
    mod_midnight = bin_mod_midnight_times(mod_midnight,tod_edges);
    
    % Fix lat thing
    for i = 1:length(lat)
        if isempty(lat{i}), lat{i} = 'unspecified'; end
    end
    
    %% Get features for soz vs not
    soz = summ.soz.chs;
    chnums = 1:length(labels);
    is_soz = ismember(chnums,soz);
    
    %% Find and remove non-intracranial
    %{
    MUST REMEMBER TO ADD THIS FOR COA
    %}
    ekg = find_non_intracranial(labels);
    
    ad = ad(~ekg,:);
    ad = nanmean(ad,1);
    
    loc = loc(~ekg,:);
    lat = lat(~ekg);
    spikes = spikes(~ekg,:); % spike rate #spikes/elec/one minute block (spikes/elec/min)
    rl = rl(~ekg,:);
    labels = labels(~ekg);
    ns = ns(~ekg,:);
    
    
    is_soz = is_soz(~ekg);
    soz_text = cell(sum(~ekg),1);
    soz_text(is_soz) = {'SOZ'};
    soz_text(~is_soz) = {'Not SOZ'};
    
    %% All rates and rl
    all_is_soz{p} = is_soz;
    all_elecs_rates{p} = nanmean(spikes,2); % mean over times
    all_elecs_rl{p} = nanmean(rl,2); % mean over times
    all_elecs_names{p} = labels;
    all_elecs_leader{p} = nansum(leader,2); % I don't think I use this
       
    %% Determine "wake" and "sleep" times
    [sleep,wake] = find_sleep_wake(ad,exc,disc);
       
    n_sleep_wake(p,1) = sum(sleep);
    n_sleep_wake(p,2) = sum(wake);
    
    %% Get wake and sleep designations for times of day
    %all_tod = 1:n_tod_bins;
    tod_sw = nan(n_tod_bins,2);
    for t = 1:n_tod_bins
        %tt = all_tod(t); % get the time of day bin im considering
        curr_bins = mod_midnight == t; % which runs match that time of day
        tod_sw(t,:) = [sum(wake(curr_bins) == 1) sum(sleep(curr_bins) == 1)]; % how many of those runs are wake and sleep
    end
    all_tod_sw(p,:,:) = tod_sw;
    
    
    
    %% wake vs sleep spike rate
    % overall spike rate (averaged across electrodes)
    mean_spikes = nanmean(spikes,1); % still spikes/elec/min
    overall_rates(p) = nanmean(spikes,'all');
    all_rates(p,:) = [nanmean(spikes(:,wake),'all') nanmean(spikes(:,sleep),'all')];
    
    % I used to do mean across electrodes first and then time
    %{
    overall_rates(p) = nanmean(spikes,'all');
    all_rates(p,:) = [nanmean(mean_spikes(wake)) nanmean(mean_spikes(sleep))];
    %}
    all_elecs_rates_sw{p} = [nanmean(spikes(:,wake),2) nanmean(spikes(:,sleep),2)];
    all_elecs_rl_sw{p} = [nanmean(rl(:,wake),2) nanmean(rl(:,sleep),2)];
    all_elecs_leader_sw{p} = [nansum(leader(:,wake),2) nansum(leader(:,sleep),2)];
    
    %% Wake vs sleep coi
    all_coi(p,:) = [nanmean(coi_global(wake)) nanmean(coi_global(sleep))];
    
    %% Wake vs sleep ns
    mean_ns = nanmean(ns,1); % node strength averaged across electrodes
    ns_sw(p,:) = [nanmean(ns(:,wake),'all') nanmean(ns(:,sleep),'all')];
    all_elecs_ns_sw{p} = [nanmean(ns(:,wake),2) nanmean(ns(:,sleep),2)];
    % I used to do mean across electrodes first and then time
    %{
    ns_sw(p,:) = [nanmean(mean_ns(wake)) nanmean(mean_ns(sleep))];
    %}
    
    %% Wake vs sleep seq info
    seq_sw(p,:) = [nanmean(seq_info(1,wake)) nanmean(seq_info(1,sleep)),...
        nanmean(seq_info(2,wake)) nanmean(seq_info(2,sleep))];

    
    %% Correlation between sleep and wake RL
    %{
    sleep_rl = nanmean(rl(:,sleep),2);
    wake_rl = nanmean(rl(:,wake),2);
    spikey = nanmean(spikes,2) > min_spikes;
    nspikey(p,1) = sum(spikey);
    nspikey(p,2) = length(spikey);
    rl_sw_corr(p) = corr(sleep_rl(spikey),wake_rl(spikey),'type','spearman','rows','pairwise');
    
    %% Rank for soz electrodes in sleep and wake
    
    spikes_for_rank = spikes;
    spikes_for_rank(isnan(spikes_for_rank)) = 0; % make nan - inf so not to screw up sort
    ranking = nan(size(spikes));
    % Loop over times
    for r = 1:size(spikes,2)
        [~,I] = sort(spikes_for_rank(:,r),'descend');
        curr_rank = 1:size(spikes,1);
        curr_rank(I) = curr_rank;
        ranking(:,r) = curr_rank;
    end
    soz_median_ranking = nanmedian(ranking(is_soz,:),1);
    soz_rank_sw(p,:) = [nanmean(soz_median_ranking(wake)),nanmean(soz_median_ranking(sleep))];
    
    %% Same ranking but RL
    rl_for_rank = rl;
    rl_for_rank(isnan(rl_for_rank)) = inf; % make nan inf so not to screw up sort
    ranking = nan(size(rl_for_rank));
    % Loop over times
    for r = 1:size(rl_for_rank,2)
        [~,I] = sort(rl_for_rank(:,r),'ascend');
        curr_rank = 1:size(rl_for_rank,1);
        curr_rank(I) = curr_rank;
        ranking(:,r) = curr_rank;
    end
    soz_median_ranking = nanmedian(ranking(is_soz,:),1);
    soz_rank_sw_rl(p,:) = [nanmean(soz_median_ranking(wake)),nanmean(soz_median_ranking(sleep))];
 %}
    
    %{
    
    %% Get spectral power for each group for locs and lats
    % Loop over loc vs lat
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
            
            %
            r_ad_ana{g}(sg,p,:) = [nanmean(spikes(ic,wake),'all') nanmean(spikes(ic,sleep),'all')];
            r_rl_ana{g}(sg,p,:) = [nanmean(rl(ic,wake),'all') nanmean(rl(ic,sleep),'all')];

            % ignoring sleep/wake diff
            rate_overall_ana{g}(sg,p) = nanmean(spikes(ic,:),'all');
            rl_overall_ana{g}(sg,p) = nanmean(rl(ic,:),'all');
            %
            
            
            % Spike rate for soz vs not
            rate_strat_ana{g}(sg,p,:) = [nanmean(spikes(ic,is_soz),'all') nanmean(spikes(ic,~is_soz),'all')];
            rl_strat_ana{g}(sg,p,:) = [nanmean(rl(ic,is_soz),'all') nanmean(rl(ic,~is_soz),'all')];
            
        end
        

    end
    %}
end

%% Remove empty pts for loc analyses
%{
missing_loc = logical(missing_loc);
r_ad_ana{1}(:,missing_loc,:) = [];
rate_overall_ana{1}(:,missing_loc) = [];
npts = npts - sum(missing_loc);
%}

%% Prep out structure
%{
out.rate_overall_ana = rate_overall_ana;
out.main = main;

out.rate_strat_ana = rate_strat_ana;
out.rl_overall_ana = rl_overall_ana;
out.rl_strat_ana = rl_strat_ana;
out.r_ad_ana = r_ad_ana;
out.r_rl_ana = r_rl_ana;
%}
out.all_rates = all_rates;
out.all_coi = all_coi;
out.ns_sw = ns_sw;
out.n_sleep_wake = n_sleep_wake;
out.names = names;
out.seq_sw = seq_sw;
%out.soz_rank_sw = soz_rank_sw;
%out.soz_rank_sw_rl = soz_rank_sw_rl;
%out.rl_sw_corr = rl_sw_corr;
out.nspikey = nspikey;
out.overall_rates = overall_rates;
out.all_is_soz = all_is_soz;
out.all_elecs_rates = all_elecs_rates;
out.all_elecs_leader = all_elecs_leader;
out.all_elecs_rl = all_elecs_rl;
out.all_elecs_names = all_elecs_names;
out.all_elecs_rl_sw = all_elecs_rl_sw;
out.all_elecs_rates_sw = all_elecs_rates_sw;
out.all_elecs_leader_sw = all_elecs_leader_sw;
out.all_elecs_ns_sw = all_elecs_ns_sw;
out.all_tod_sw = all_tod_sw;
out.tod_edges = tod_edges;

%% (No sleep) How does spike rate and timing vary across locations
%{
f1 = figure;
set(gcf,'position',[10 271 1260 526])
tiledlayout(2,3)

% spike rate by location
nexttile
curr_rate = rate_overall_ana{1}; % location
plot_paired_data(curr_rate,main_locs,'Spike/elec/min','paired')

% spike rate by soz vs not
nexttile
plot_paired_data(rate_soz',main_soz,'Spike/elec/min','paired')

% Is spike rate higher in SOZ within each anatomical region
nexttile
curr_rate = rate_strat_ana{1};
interaction_plot_and_stats(curr_rate*1e3,main_locs,...
    'Spikes/elec/min',{'SOZ','Not SOZ'},1);

% Spike rl by location
nexttile
curr_rl = rl_overall_ana{1}; % location
plot_paired_data(curr_rl*1e3,main_locs,'Spike latency (ms)','paired')

% spike rl by SOZ
nexttile
plot_paired_data(rl_soz'*1e3,main_soz,'Spike latency (ms)','paired')


% Is rl lower in SOZ within each anatomical region
nexttile
curr_rl = rl_strat_ana{1};
interaction_plot_and_stats(curr_rl*1e3,main_locs,...
    'Spike latency (ms)',{'SOZ','Not SOZ'},1);
print(f1,[out_folder,'no_sleep'],'-dpng')

%% How do spikes vary with sleep
f2 = figure;
set(gcf,'position',[10 10 800 1000])
tiledlayout(3,2,'tilespacing','tight','padding','tight')

% ROC
nexttile
plot(roc(:,1),roc(:,2),'k','linewidth',2)
hold on
plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','northwest')
set(gca,'fontsize',15)

% Overall spike rate sleep vs wake
nexttile
plot_paired_data(all_rates',{'Wake','Sleep'},'Spike/elec/min','paired')

% COI sleep vs wake
nexttile
plot_paired_data(all_coi',{'Wake','Sleep'},'Spike COI','paired')

% spike rate consistency wake vs sleep
nexttile
plot_paired_data(all_src',{'Wake','Sleep'},'Spike rate consistency','paired')

% spike timing consistency wake vs sleep
nexttile
plot_paired_data(all_stc',{'Wake','Sleep'},'Spike timing consistency','paired')

% average ns wake vs sleep
nexttile
plot_paired_data(ns_sw',{'Wake','Sleep'},'Average node strength','paired')
print(f2,[out_folder,'sleep_fig'],'-dpng')



%% Interaction between sleep and location
f3=figure;
set(gcf,'position',[10 10 600 500])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

% Is sleep-related increase in spike rate higher for SOZ?
nexttile
interaction_plot_and_stats(rate_sw_soz,main_soz,'Spike/elec/min',{'Wake','Sleep'},0);

nexttile
plot_and_stats_change(rate_sw_soz,main_soz,{'Spike rate change','in sleep'},'paired')

% Is sleep-related increase in spike rate higher for different anatomical
% locations?
nexttile
interaction_plot_and_stats(r_ad_ana{1},main_locs,'Spike/elec/min',{'Wake','Sleep'},0);

nexttile
plot_and_stats_change(r_ad_ana{1},main_locs,{'Spike rate change','in sleep'},'paired')
print(f3,[out_folder,'sleep_loc_interaction'],'-dpng')
%}

end