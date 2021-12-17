function sleep_figure2

%% Parameters
plot_type = 'scatter';
nblocks = 6;
myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];



locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];


%% Load out file and get roc stuff
out = load([out_folder,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

figure
set(gcf,'position',[100 100 900 700])
tiledlayout(4,2,'tilespacing','tight','padding','tight')

%% Post-ictal surge
all_pts_spikes_bins = sz_out.all_pts_spikes_bins;
all_pts_sleep_bins = sz_out.all_pts_sleep_bins;
surround_hours = sz_out.surround_hours;

ax1 = nexttile;
npts = size(all_pts_spikes_bins,1);
median_spikes = nanmedian(all_pts_spikes_bins,1);
mean_sleep = nanmean(all_pts_sleep_bins,1)*100;
iqr_spikes = prctile(all_pts_spikes_bins,[25,75],1);
times = linspace(-surround_hours,surround_hours,length(median_spikes));
%ylabels = ["Spikes/elecs/min","% Asleep"];
shaded_error_bars(times,median_spikes,iqr_spikes,myColours(1,:));
hold on
plot([0 0],ylim,'k--','linewidth',2)

%% Pre vs post ictal
nexttile([2 1])
nbins = size(all_pts_spikes_bins,2);
pre = 1:nbins/2;
post = nbins/2+1:nbins;
pre_post = [nanmean(all_pts_spikes_bins(:,pre),2),nanmean(all_pts_spikes_bins(:,post),2)];
plot_paired_data(pre_post',{'pre-ictal','post-ictal'},'Spikes/elec/min','paired',plot_type)

%% Percent detected asleep
nexttile
plot(times,mean_sleep,'color',myColours(2,:),'linewidth',2)
hold on
plot([0 0],ylim,'k--','linewidth',2)

%% Temporal vs extratemporal
nexttile([2 1])
loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
rate_diff = (pre_post(:,2)-pre_post(:,1));
p = ranksum(rate_diff(temporal),rate_diff(extra));
plot(1+randn(sum(temporal),1)*0.05,rate_diff(temporal),'o','linewidth',2,'color',myColours(1,:))
hold on
plot([0.7 1.3],[nanmedian(rate_diff(temporal)) nanmedian(rate_diff(temporal))],...
    'linewidth',2,'color',myColours(1,:))
plot(2+randn(sum(extra),1)*0.05,rate_diff(extra),'o','linewidth',2,'color',myColours(2,:))
plot([1.7 2.3],[nanmedian(rate_diff(extra)) nanmedian(rate_diff(extra))],...
    'linewidth',2,'color',myColours(2,:))
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})
ylabel('Sleep-wake spikes/elec/min')
title('Sleep-wake spike rate difference')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)

end