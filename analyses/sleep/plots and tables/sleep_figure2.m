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
set(gcf,'position',[100 100 1400 700])
tiledlayout(4,3,'tilespacing','compact','padding','compact')

%% Seizure timing PSD
nexttile([2 1])
periods = sz_circ_out.periods;
iqr_psd = sz_circ_out.iqr_psd;
all_psd = sz_circ_out.all_psd;
median_psd = nanmedian(all_psd,1);
mp = shaded_error_bars(periods,median_psd,iqr_psd,[0 0 0]);
set(gca,'fontsize',15)
xlabel('Period (hours)')
ylabel('Power index')
title('Seizure time periodogram')

%% Percent asleep
nexttile([2 1])
pre_wake = sz_circ_out.pre_wake;
n_sleep_wake = bin_out.n_sleep_wake;
perc_sz_asleep = cellfun(@(x) prc_asleep(x),pre_wake);
perc_all_asleep = 100*n_sleep_wake(:,1)./sum(n_sleep_wake,2);
minp = min([perc_sz_asleep;perc_all_asleep]);
maxp = max([perc_sz_asleep;perc_all_asleep]);
plot(perc_all_asleep,perc_sz_asleep,'ko','linewidth',2)
hold on
%plot([0 100],[0 100],'k--','linewidth',2)
xlim([0 100])
ylim([0 100])
plot([nanmedian(perc_all_asleep) nanmedian(perc_all_asleep)],ylim,'--','color','k','linewidth',2)
text(nanmedian(perc_all_asleep),91,sprintf('\\leftarrow Median %1.1f%% asleep',nanmedian(perc_all_asleep)),...
    'fontsize',15,'color','k')

xlabel('Total time asleep (%)')
ylabel('Seizures arising from sleep (%)')
set(gca,'fontsize',15)
title('Percentage of all times and seizures from sleep')

nexttile([2 1])
loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
p = ranksum(perc_sz_asleep(temporal),perc_sz_asleep(extra));
plot(1+randn(sum(temporal),1)*0.05,perc_sz_asleep(temporal),'o','linewidth',2,'color',myColours(1,:))
hold on
plot([0.7 1.3],[nanmedian(perc_sz_asleep(temporal)) nanmedian(perc_sz_asleep(temporal))],...
    'linewidth',2,'color',myColours(1,:))
plot(2+randn(sum(extra),1)*0.05,perc_sz_asleep(extra),'o','linewidth',2,'color',myColours(2,:))
plot([1.7 2.3],[nanmedian(perc_sz_asleep(extra)) nanmedian(perc_sz_asleep(extra))],...
    'linewidth',2,'color',myColours(2,:))
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})
ylabel('Seizures arising from sleep (%)')
title('Seizure state-dependence by localization')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)


%% Post-ictal surge
all_pts_spikes_bins = sz_out.all_pts_spikes_bins;
all_pts_sleep_bins = sz_out.all_pts_sleep_bins;
surround_hours = sz_out.surround_hours;

ax1 = nexttile([1 1]);
npts = size(all_pts_spikes_bins,1);
median_spikes = nanmedian(all_pts_spikes_bins,1);
mean_sleep = nanmean(all_pts_sleep_bins,1)*100;
iqr_spikes = prctile(all_pts_spikes_bins,[25,75],1);
times = linspace(-surround_hours,surround_hours,length(median_spikes));
%ylabels = ["Spikes/elecs/min","% Asleep"];
shaded_error_bars(times,median_spikes,iqr_spikes,myColours(1,:));
hold on
plot([0 0],ylim,'k--','linewidth',3)
ylabel('Spikes/elecs/min')
title('Peri-ictal spike rates')
set(gca,'fontsize',15)
xlabel('Hours surrounding seizure')

%% Pre vs post ictal
nexttile([2 1])
nbins = size(all_pts_spikes_bins,2);
pre = 1:nbins/2;
post = nbins/2+1:nbins;
pre_post = [nanmean(all_pts_spikes_bins(:,pre),2),nanmean(all_pts_spikes_bins(:,post),2)];
plot_paired_data(pre_post',{'pre-ictal state','post-ictal state','post-ictally'},'Spikes/elec/min','paired',plot_type)
title('Pre- vs post-ictal spike rates')


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
ylabel('Post-pre-ictal spikes/elec/min')
title('Post-pre-ictal spike rate difference by localization')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)

%% Percent detected asleep
nexttile([1 1])
plot(times,mean_sleep,'color',myColours(2,:),'linewidth',2)
hold on
plot([0 0],ylim,'k--','linewidth',3)
set(gca,'fontsize',15)
ylabel('% Classified asleep')
title('Peri-ictal sleep classification')
xlabel('Hours surrounding seizure')

%% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.32 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0.64 0.91 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.42 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.19 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0.32 0.42 0.1 0.1],'String','F','fontsize',25,'linestyle','none')
annotation('textbox',[0.64 0.42 0.1 0.1],'String','G','fontsize',25,'linestyle','none')



print([out_folder,'Fig2'],'-dpng')

end

function prc = prc_asleep(x)

prc = 100 * sum(x==0)/(sum(x==1)+sum(x==0));

end