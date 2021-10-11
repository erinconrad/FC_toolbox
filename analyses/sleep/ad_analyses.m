function ad_analyses(summ)

%% To do
%{
- exclude sz times
- get soz to do that analysis
%}

%% Parameters
main_locs = {'mesial temporal','temporal neocortical','other cortex','white matter'};
main_lats = {'Left','Right'};
main{1} = main_locs;
main{2} = main_lats;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
summ_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Alpha delta ratio validation
swdes = sw_ad_erin_designations;
npts_val = length(swdes);
ad_norm = nan(npts_val,2); %1 = sleep, 2 = wake
all_wake = [];
all_sleep = [];
for j = 1:npts_val
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
[roc,auc] = calculate_roc(all_sleep,all_wake,1e3);

%% Main analyses
npts = length(summ);
r_ad_spikes = nan(npts,1);
skip_pts = [];

r_ad_ana = cell(2,1);
for i = 1:length(r_ad_ana)
    r_ad_ana{i} = nan(length(main{i}),npts);
end

for p = 1:npts
    
    %% Get main things
    loc = summ(p).ana_loc;
    lat = summ(p).ana_lat;
    spikes = summ(p).spikes;
    ad = summ(p).ad;
    ad = nanmean(ad,1);
    
    
    % Skip if all empty
    if sum(cellfun(@(x) isempty(x),loc)) == length(loc) 
        fprintf('\nskipping pt %d\n',p);
        skip_pts = [skip_pts;p];
        continue
    end
    
    
    %% Correlation between spike rate and ad
    % overall spike rate (averaged across electrodes)
    mean_spikes = nanmean(spikes,1);
    r_ad_spikes(p) = corr(mean_spikes',ad','rows','pairwise');
    
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
            rate_subgroup = nanmean(spikes(ic,:),1);
                    
            % Get the spike rate/ad correlation for that region
            r_ad_ana{g}(sg,p) = corr(rate_subgroup',ad','rows','pairwise');

        end
        

    end
end

%% Remove empty pts
r_ad_spikes(skip_pts) = [];
for i = 1:length(r_ad_ana)
    r_ad_ana{i}(:,skip_pts) = [];
end
npts = npts - length(skip_pts);


%% initialize figure
figure
set(gcf,'position',[100 100 800 600])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

%% A/D ratio in sleep vs wake
nexttile
for p = 1:size(ad_norm,1)
    if any(isnan(ad_norm(p,:))), continue; end
    plot([1 2],[ad_norm(p,1),ad_norm(p,2)],'k-');
    hold on
end
xticks([1 2])
xticklabels({'Sleep','Awake'})
ylabel({'Alpha/delta ratio','(normalized)'})
xlim([0.5 2.5])
set(gca,'fontsize',15)
[~,pval] = ttest(ad_norm(:,1),ad_norm(:,2));

%% ROC curve
nexttile
plot(roc(:,1),roc(:,2),'k','linewidth',2)
hold on
plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','northwest')
set(gca,'fontsize',15)

%% All corr coeffs between AD ratio and spikes
nexttile
plot(r_ad_spikes,'ko','linewidth',2)
hold on
xlim([0 npts+1])
plot(xlim,[0 0],'k--')
ylim([-1 1])
set(gca,'fontsize',15)
xticklabels([])
xlabel('Patient')
ylabel('AD-spike rate correlation')

% Fisher transform rs
z_ad_spikes = atanh(r_ad_spikes); 
% two sided ttest of the fisher transformed rs
[~,pval,~,stats] = ttest(z_ad_spikes);
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('%s',get_p_text(pval)),...
    'fontsize',15,'verticalalignment','top')

%% AD/spike corr by anatomical locations
nexttile
loc_corr = r_ad_ana{1}; % loc
median_corr = nanmedian(loc_corr,2);
iqr_corr = [prctile(loc_corr,25,2),prctile(loc_corr,75,2)];
errorbar(1:length(main_locs),median_corr ,...
    median_corr -iqr_corr(:,1),...
    iqr_corr(:,2)-median_corr ,'o','markersize',10,...
    'linewidth',2,'color','k');
hold on
ylim([-1 1])
xlim([0 length(median_corr)+1])
plot(xlim,[0 0],'k--')
xticks(1:length(median_corr))
xticklabels(main{1})
%xtickangle(30)
ylabel('AD-spike rate correlation')
set(gca,'fontsize',15)

% Do stats
[p,post_hoc_p,which_groups] = non_para_circ_stats(loc_corr);
if p > 0.05
    pairs_to_plot = [];
else
    pairs_to_plot = which_groups(post_hoc_p < 0.05/size(which_groups,1),:);
    post_hoc_p_to_plot = post_hoc_p(post_hoc_p < 0.05/size(which_groups,1));
end
yl=ylim;
heights = get_heights([-1 0.3],pairs_to_plot);
%ylim([yl(1) heights(end,2)]);

if p > 0.05
    plot([1 length(avg_over_pts)],...
        [heights(size(heights,1)-1,1) heights(size(heights,1)-1,1)],'k',...
        'linewidth',2);
    text(mean([1 length(avg_over_pts)]),heights(size(heights,1)-1,2),...
        'ns','fontsize',15,'horizontalalignment','center')
else
    for k = 1:size(pairs_to_plot,1)
        plot([pairs_to_plot(k,1)+0.1 pairs_to_plot(k,2)-0.1],[heights(k,1) heights(k,1)],'k-',...
            'linewidth',2)
        hold on
        text(mean(pairs_to_plot(k,:)),heights(k,2),...
            get_asterisks(post_hoc_p_to_plot(k),size(which_groups,1)),...
            'fontsize',15,'horizontalalignment','center')
    end
end

print([out_folder,'ad_analyses'],'-dpng')
close all
end