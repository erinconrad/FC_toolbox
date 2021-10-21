function binary_ad_analyses

%% To do
%{
- exclude sz times

- is overlap between most frequent or earliest spiker and soz higher in
sleep


%}

%% Questions
%{
Fig 1 (no sleep/wake dependence) - 7 plots
- spike rate by location and laterality and SOZ (3 plots)
- rl by location and laterality and SOZ (3 plots)
- spike rate correlates with rl?

Fig 2 (sleep wake, ignoring antomy) -  5 plots
- ROC curve for alpha delta ratio
- overall spike rate sleep vs wake
- COI sleep vs wake
- spike rate consistency wake vs sleep
- spike timing consistency wake vs sleep

Fig 3 (sleep wake/anatomy interaction)
- rate sleep vs wake by location and laterality and SOZ
- rl sleep vs wake by location and laterality and SOZ

- How does overall spike rate change with sleep
- Does the correlation between sleep and spike rate depend on anatomical
location?
- How does spike spread change with sleep?
- Does spike timing correlate with spike rate? Is the spikiest channel also
the one that spikes first?
- Is the order of spike rate more consistent in wake or sleep?
- Is the order of spike timing more consistent in wake or sleep?
     - For these, I would probably need to designate wake and sleep
%}

%% Parameters
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

%% Alpha delta ratio validation
fprintf('\nDoing alpha delta ratio validation\n');

swdes = sw_ad_erin_designations;
npts_val = length(swdes);
ad_norm = nan(npts_val,2); %1 = sleep, 2 = wake
all_wake = [];
all_sleep = [];
all_rate_rl_corr = [];
for j = 1:npts_val
    if isempty(swdes(j).sw), continue; end
    sleep_ad = swdes(j).sw.sleep;
    wake_ad = swdes(j).sw.wake;
    ad_val = swdes(j).ad;

    sleep_norm = (sleep_ad-nanmedian(ad_val))./iqr(ad_val);
    wake_norm = (wake_ad-nanmedian(ad_val))./iqr(ad_val);
    ad_norm(j,:) = [nanmean(sleep_norm),nanmean(wake_norm)];
    all_wake = [all_wake;wake_norm];
    all_sleep = [all_sleep;sleep_norm];
end

% Calculate roc
[roc,auc,disc] = calculate_roc(all_sleep,all_wake,1e3);


fprintf('\nGot alpha delta ratio validation\n');

%% Main analyses
missing_loc = [];
all_rates = [];
all_coi = [];
all_src = [];
all_stc = [];

r_ad_ana = cell(2,1);
for i = 1:length(r_ad_ana)
    r_ad_ana{i} = nan(length(main{i}),npts,2);
end

r_rl_ana = cell(2,1);
for i = 1:length(r_rl_ana)
    r_rl_ana{i} = nan(length(main{i}),npts,2);
end

rate_overall_ana = cell(2,1);
for i = 1:length(rate_overall_ana)
    rate_overall_ana{i} = nan(length(main{i}),npts);
end

rl_overall_ana = cell(2,1);
for i = 1:length(rl_overall_ana)
    rl_overall_ana{i} = nan(length(main{i}),npts);
end


rate_soz = nan(npts,2);
rl_soz = nan(npts,2);
rate_sw_soz = nan(2,npts,2); % 2 for soz/not, then 2 for wake/sleep
rl_sw_soz = nan(2,npts,2);% 2 for soz/not, then 2 for wake/sleep
overlap_spikiest = nan(npts,2); % wake/sleep
overlap_earliest = nan(npts,2); % wake/sleep
null_ps = nan(npts,1);
ns_sw = nan(npts,2);
ge_sw = nan(npts,2);

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
    ge = summ.ge;
    
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
    spikes = spikes(~ekg,:);
    rl = rl(~ekg,:);
    labels = labels(~ekg);
    ns = ns(~ekg,:);
    
    
    is_soz = is_soz(~ekg);
       
    %% Determine "wake" and "sleep" times
    % normalized ad
    ad_norm = (ad - nanmedian(ad))./iqr(ad);
    wake = ad_norm > disc;
    sleep = ad_norm <= disc;
       
    
    %% wake vs sleep spike rate
    % overall spike rate (averaged across electrodes)
    mean_spikes = nanmean(spikes,1);
    all_rates = [all_rates;nanmean(mean_spikes(wake)) nanmean(mean_spikes(sleep))];
    
    %% Wake vs sleep coi
    all_coi = [all_coi;nanmean(coi_global(wake)) nanmean(coi_global(sleep))];
    
    %% Wake vs sleep ns
    mean_ns = nanmean(ns,1); % node strength averaged across electrodes
    ns_sw = [ns_sw;nanmean(mean_ns(wake)) nanmean(mean_ns(sleep))];
    ge_sw = [ge_sw;nanmean(ge(wake)) nanmean(ge(sleep))];

    
    %% SRC - spike rate consistency
    % Spikes in wake and sleep
    wake_spikes = spikes(:,wake);
    sleep_spikes = spikes(:,sleep);
    
    % Mean vector of spike rates across electrodes
    mean_wake_spikes = nanmean(wake_spikes,2);
    mean_sleep_spikes = nanmean(sleep_spikes,2);
    
    % SRC - Spearman correlation of spikes rates with the mean
    src_wake = nanmean(corr(mean_wake_spikes,wake_spikes,'type','spearman',...
        'rows','pairwise'));
    src_sleep = nanmean(corr(mean_sleep_spikes,sleep_spikes,'type','spearman',...
        'rows','pairwise'));
    all_src = [all_src;src_wake src_sleep];
    
    
    %% Rate-rl correlation
    avg_rate = nanmean(spikes,2);
    avg_rl = nanmean(rl,2);
    rate_rl_corr = corr(avg_rate,avg_rl,'type','spearman','rows','pairwise');
    all_rate_rl_corr = [all_rate_rl_corr;rate_rl_corr];
    
    %% SRC - spike timing consistency
    % Spikes in wake and sleep
    wake_rl = rl(:,wake);
    sleep_rl = rl(:,sleep);
    
    % Mean vector of spike timing (rl) across electrodes
    mean_wake_rl = nanmean(wake_rl,2);
    mean_sleep_rl = nanmean(sleep_rl,2);
    
    % SRC - Spearman correlation of spikes timing with the mean
    stc_wake = nanmean(corr(mean_wake_rl,wake_rl,'type','spearman',...
        'rows','pairwise'));
    stc_sleep = nanmean(corr(mean_sleep_rl,sleep_rl,'type','spearman',...
        'rows','pairwise'));
    all_stc = [all_stc;stc_wake stc_sleep];
    
    %% Spikiest and earliest in sleep and wake    
    [~,spikiest_sw(1)] = max(mean_wake_spikes);
    [~,spikiest_sw(2)] = max(mean_sleep_spikes);
    
    [~,earliest_sw(1)] = min(mean_wake_rl);
    [~,earliest_sw(2)] = min(mean_sleep_rl);
    
    % Do these overlap with SOZ
    overlap_spikiest(p,:) = [is_soz(spikiest_sw(1)) is_soz(spikiest_sw(2))];
    overlap_earliest(p,:) = [is_soz(earliest_sw(1)) is_soz(earliest_sw(2))];
    null_ps(p) = sum(is_soz)/length(is_soz);
    
    %% SOZ analysis
    % average over all times, all soz (and separately for all non soz)
    rate_soz(p,:) = [nanmean(spikes(is_soz,:),'all') nanmean(spikes(~is_soz,:),'all')];
    rl_soz(p,:) = [nanmean(rl(is_soz,:),'all') nanmean(rl(~is_soz,:),'all')];
    
    % Now sleep vs wake
    for sz = 1:2
        % first SOZ
        if sz == 1
            rate_sw_soz(sz,p,:) = [nanmean(spikes(is_soz,wake),'all') nanmean(spikes(is_soz,sleep),'all')];
            rl_sw_soz(sz,p,:) = [nanmean(rl(is_soz,wake),'all') nanmean(rl(is_soz,sleep),'all')];
        else
            % not soz
            rate_sw_soz(sz,p,:) = [nanmean(spikes(~is_soz,wake),'all') nanmean(spikes(~is_soz,sleep),'all')];
            rl_sw_soz(sz,p,:) = [nanmean(rl(~is_soz,wake),'all') nanmean(rl(~is_soz,sleep),'all')];
        end
    end
    
    %% Skip subsequent loc analyses if missing
    if sum(cellfun(@(x) isempty(x),loc)) == length(loc) 
        missing_loc = [missing_loc;p];
        continue
    end
    
 
    %% Get spectral power for each group for locs and lats
    % Loop over loc vs lat
    for g = 1:2
        if g == 1
            group = loc;
        elseif g == 2
            group = lat;
        end
        
        % Get the rates corresponding to the subgroups
        % (can probably do this without a for loop)
       
        for sg = 1:length(main{g})
            ic = ismember(group,main{g}(sg));
            
            %
            % Get the spike rate/ad correlation for that region
            r_ad_ana{g}(sg,p,:) = [nanmean(spikes(ic,wake),'all') nanmean(spikes(ic,sleep),'all')];
            r_rl_ana{g}(sg,p,:) = [nanmean(rl(ic,wake),'all') nanmean(rl(ic,sleep),'all')];

            % ignoring sleep/wake diff
            rate_overall_ana{g}(sg,p) = nanmean(spikes(ic,:),'all');
            rl_overall_ana{g}(sg,p) = nanmean(rl(ic,:),'all');
            %}
            
            %{
            rate_subgroup = nanmean(spikes(ic,:),1);
            rl_subgroup = nanmean(rl(ic,:),1);
                    
            % Get the spike rate/ad correlation for that region
            r_ad_ana{g}(sg,p,:) = [nanmean(rate_subgroup(wake)) nanmean(rate_subgroup(sleep))];
            r_rl_ana{g}(sg,p,:) = [nanmean(rl_subgroup(wake)) nanmean(rl_subgroup(sleep))];

            % ignoring sleep/wake diff
            rate_overall_ana{g}(sg,p) = nanmean(rate_subgroup);
            rl_overall_ana{g}(sg,p) = nanmean(rl_subgroup);
            %}
        end
        

    end
end

%% Remove empty pts for loc analyses
for i = 1:length(r_ad_ana)
    r_ad_ana{i}(:,missing_loc,:) = [];
    rate_overall_ana{i}(:,missing_loc) = [];
end
npts = npts - length(missing_loc);


%% Figure 1 - ignoring sleep, comparison across locations
%{
Fig 1 (no sleep/wake dependence) - 7 plots
- spike rate by location and laterality and SOZ (3 plots)
- rl by location and laterality and SOZ (3 plots)
- spike rate correlates with rl?
%}
f1 = figure;
set(f1,'position',[100 100 800 600])
tiledlayout(3,3,'tilespacing','compact','padding','compact')

% rate overall location (how does spike rate vary by location)
nexttile
curr_rate = rate_overall_ana{1}; % location
plot_paired_data(curr_rate,main_locs,'Spike/elec/???')

% rloverall location
nexttile
curr_rl = rl_overall_ana{1}; % location
plot_paired_data(curr_rl*1e3,main_locs,'Spike latency (ms)')

%rl-rate correlation
nexttile([3 1])
plot(all_rate_rl_corr,'o','linewidth',2)
xlim([0 length(all_rate_rl_corr)+1])
hold on
plot(xlim,[0 0],'k--')
ylabel('Spike rate-latency correlation')
z = fisher_transform(all_rate_rl_corr);
[~,p] = ttest(z);
simple_p_stats(ylim,p,[1 length(z)])
set(gca,'fontsize',15)


% rate overall lateralization
nexttile
curr_rate = rate_overall_ana{2}; % lateralization
plot_paired_data(curr_rate,main_lats,'Spike/elec/???')

% rloverall lateralization
nexttile
curr_rl = rl_overall_ana{2}; % lateralization
plot_paired_data(curr_rl*1e3,main_lats,'Spike latency (s)')

% Insert SOZ-rate comparison here
nexttile
plot_paired_data(rate_soz',main_soz,'Spike/elec/???')

% Insert SOZ-rl comparison here
nexttile
plot_paired_data(rl_soz'*1e3,main_soz,'Spike latency (s)')

print(f1,[out_folder,'ad_fig1'],'-dpng')


%% Figure 2
%{
Fig 2 (sleep wake, ignoring antomy)
- ROC curve for alpha delta ratio
- overall spike rate sleep vs wake
- COI sleep vs wake
- spike rate consistency wake vs sleep
- spike timing consistency wake vs sleep
%}
f2 = figure;
set(gcf,'position',[100 100 900 700])
tiledlayout(3,2,'tilespacing','compact','padding','compact')

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
plot_paired_data(all_rates',{'Wake','Sleep'},'Spike/elec/???')

% COI sleep vs wake
nexttile
plot_paired_data(all_coi',{'Wake','Sleep'},'Spike COI')

% spike rate consistency wake vs sleep
nexttile
plot_paired_data(all_src',{'Wake','Sleep'},'Spike rate consistency')

% spike timing consistency wake vs sleep
nexttile
plot_paired_data(all_stc',{'Wake','Sleep'},'Spike timing consistency')

% GE wake vs sleep
nexttile
plot_paired_data(ge_sw',{'Wake','Sleep'},'Global efficiency')
print(f2,[out_folder,'ad_fig2'],'-dpng')

%% Figure 3 - sleep/location interaction
%{
Fig 3 (sleep wake/anatomy interaction)
- rate sleep vs wake by location and laterality and SOZ
- rl sleep vs wake by location and laterality and SOZ
%}
f3 = figure;
set(gcf,'position',[100 100 800 800])
tiledlayout(3,2,'tilespacing','compact','padding','compact')

% Rate sleep vs wake by location
nexttile
interaction_plot_and_stats(r_ad_ana{1},main_locs,'Spike/elec/???',{'Wake','Sleep'},0);

% Rl sleep vs wake by location -  this is cool. It appears that the latency
% doesn't go up much with sleep for mesial temporal, but it does for other
% locations. That suggests that in sleep, spikes start in the same place
% they always started (the mesial temporal lobe) and then spread to other
% areas more.
nexttile
interaction_plot_and_stats(r_rl_ana{1}*1e3,main_locs,'Spike latency (ms)',{'Wake','Sleep'},1);

% Rate sleep vs wake by laterality
%{
nexttile
interaction_plot_and_stats(r_ad_ana{2},main_lats,'Spike/elec/???',{'Wake','Sleep'},0);

% Rl sleep vs wake by laterality
nexttile
interaction_plot_and_stats(r_rl_ana{2}*1e3,main_lats,'Spike latency (ms)',{'Wake','Sleep'},1);
%}

% Rate sleep vs wake by SOZ
nexttile
interaction_plot_and_stats(rate_sw_soz,main_soz,'Spike/elec/???',{'Wake','Sleep'},0);

% Rl sleep vs wake by SOZ
nexttile
interaction_plot_and_stats(rl_sw_soz*1e3,main_soz,'Spike latency (ms)',{'Wake','Sleep'},1);

% Overlap with highest spiker and SOZ in wake vs sleep
nexttile
plot_overlap(overlap_spikiest,{'Wake','Sleep'},...
    {'How often is spikiest','electrode in SOZ'},1e4,null_ps);

% Overlap with earliest spiker and SOZ in wake vs sleep
nexttile
plot_overlap(overlap_earliest,{'Wake','Sleep'},...
    {'How often is earliest','spiking electrode in SOZ'},1e4,null_ps);


print(f3,[out_folder,'ad_fig3'],'-dpng')
%{
if 0
    
%% Figure 2 - changes with sleep
figure
set(gcf,'position',[100 100 800 600])
tiledlayout(3,2,'tilespacing','compact','padding','compact')

%% ROC curve
%{
nexttile
plot(roc(:,1),roc(:,2),'k','linewidth',2)
hold on
plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','northwest')
set(gca,'fontsize',15)
%}

%% Spikes in sleep vs wake
nexttile
for p = 1:size(all_rates,1)
    if any(isnan(all_rates(p,:))), continue; end
    plot([1 2],[all_rates(p,1),all_rates(p,2)],'k-');
    hold on
end
xticks([1 2])
xticklabels({'Awake','Sleep'})
ylabel({'Spikes/elec/???'})
xlim([0.5 2.5])
set(gca,'fontsize',15)
[~,pval] = ttest(all_rates(:,1),all_rates(:,2));
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
plot([1 2],[ybar ybar],'k-');
text(1.5,ytext,get_asterisks(pval,1));

%% COI in sleep vs wake
nexttile
for p = 1:size(all_coi,1)
    if any(isnan(all_coi(p,:))), continue; end
    plot([1 2],[all_coi(p,1),all_coi(p,2)],'k-');
    hold on
end
xticks([1 2])
xticklabels({'Awake','Sleep'})
ylabel({'Spike COI'})
xlim([0.5 2.5])
set(gca,'fontsize',15)
[~,pval] = ttest(all_coi(:,1),all_coi(:,2));
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
plot([1 2],[ybar ybar],'k-');
text(1.5,ytext,get_asterisks(pval,1));

%% SRC in sleep vs wake
nexttile
for p = 1:size(all_src,1)
    if any(isnan(all_src(p,:))), continue; end
    plot([1 2],[all_src(p,1),all_src(p,2)],'k-');
    hold on
end
xticks([1 2])
xticklabels({'Awake','Sleep'})
ylabel({'Spike rate consistency'})
xlim([0.5 2.5])
set(gca,'fontsize',15)
[~,pval] = ttest(all_src(:,1),all_src(:,2));
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
plot([1 2],[ybar ybar],'k-');
text(1.5,ytext,get_asterisks(pval,1));

%% STC in sleep vs wake
nexttile
for p = 1:size(all_stc,1)
    if any(isnan(all_stc(p,:))), continue; end
    plot([1 2],[all_stc(p,1),all_stc(p,2)],'k-');
    hold on
end
xticks([1 2])
xticklabels({'Awake','Sleep'})
ylabel({'Spike timing consistency'})
xlim([0.5 2.5])
set(gca,'fontsize',15)
[~,pval] = ttest(all_stc(:,1),all_stc(:,2));
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
newyl = yl(1) + 1.15*(yl(2)-yl(1));
plot([1 2],[ybar ybar],'k-');
ylim([yl(1) newyl]);
text(1.5,ytext,get_asterisks(pval,1));

%% Spike rate in sleep vs wake by anatomical locations
nexttile
loc_rates = r_ad_ana{1}; % loc

% loop over anatomical locations
for i = 1:size(loc_rates,1)
    errorbar([i-0.2 i+0.2],squeeze(nanmean(loc_rates(i,:,:),2)),...
        squeeze(nanstd(loc_rates(i,:,:),[],2)),'o')
    hold on
end
xticks(1:size(loc_rates,1))
xticklabels(main_locs)

yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
newyl = yl(1) + 1.15*(yl(2)-yl(1));
for i = 1:size(loc_rates,1)
    [~,pval] = ttest(loc_rates(i,:,1),loc_rates(i,:,2));
    plot([i-0.2 i+0.2],[ybar ybar],'k-');
    text(i,ytext,get_asterisks(pval,size(loc_rates,1)));
end
ylim([yl(1) newyl])

%% Spike rate in sleep vs wake by anatomical locations
nexttile
loc_rl = r_rl_ana{1}; % loc

% loop over anatomical locations
for i = 1:size(loc_rl,1)
    errorbar([i-0.2 i+0.2],squeeze(nanmean(loc_rl(i,:,:),2)),...
        squeeze(nanstd(loc_rl(i,:,:),[],2)),'o')
    hold on
end
xticks(1:size(loc_rl,1))
xticklabels(main_locs)

yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
newyl = yl(1) + 1.15*(yl(2)-yl(1));
for i = 1:size(loc_rates,1)
    [~,pval] = ttest(loc_rl(i,:,1),loc_rl(i,:,2));
    plot([i-0.2 i+0.2],[ybar ybar],'k-');
    text(i,ytext,get_asterisks(pval,size(loc_rl,1)));
end
ylim([yl(1) newyl])




end
    %}





close all
end