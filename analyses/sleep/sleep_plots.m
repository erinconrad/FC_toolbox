function sleep_plots(out)

%% Parameters
nblocks = 6;
baseline = 1:nblocks; % one hour
preictal = 31:31+nblocks-1;
postictal = 37:37+nblocks-1;

%% Unpack substructures
unpack_any_struct(out);

%% Fig 1 - Circadian analysis
figure
set(gcf,'position',[100 100 1400 500])
tiledlayout(1,3,'tilespacing','tight','padding','tight')

% A - PSD
median_psd = circ_out.median_psd;
iqr_psd = circ_out.iqr_psd;
periods = circ_out.periods;
low_period = circ_out.low_period;
circ_P = circ_out.circ_P;
all_circ_P = circ_out.all_circ_P;
main_locs = circ_out.main_locs;
all_locs = circ_out.all_locs;

nexttile
shaded_error_bars(periods,median_psd,iqr_psd,[0 0 0]);
xlim([0 100])
xlabel('Period (hours)')
ylabel({'Spike rate', 'normalized power spectrum'});
set(gca,'fontsize',15)
title('Spike rate power spectral density')

% Compare cyclical power across electrode localizations
nexttile
plot_paired_data(circ_P{1},main_locs,'Relative circadian power','paired','violin')
title({'Relative circadian power','by spike location'})

% Compare cyclical power between patients with different soz localizations
nexttile
mt = cellfun(@(x) strcmp(x,'temporal'), all_locs);
circ_mt = all_circ_P(mt);
circ_other = all_circ_P(~mt);
circ_mt = [circ_mt;nan(length(mt)-length(circ_mt),1)]; % pad both with nans
circ_other = [circ_other;nan(length(mt)-length(circ_other),1)];
plot_paired_data(([circ_mt,circ_other])',{'Temporal','Other'},'Relative circadian power',...
    'unpaired','violin')
title({'Relative circadian power','by SOZ localization'});
print([out_folder,'Fig1'],'-dpng')
close(gcf)

%% Figure 2 - What happens to spikes with sleep
figure
set(gcf,'position',[10 10 800 1000])
tiledlayout(3,2,'tilespacing','tight','padding','tight')

% A: ROC curve
roc = roc_out.roc;
auc = roc_out.auc;
disc_I = roc_out.disc_I;

nexttile([1 2])
plot(roc(:,1),roc(:,2),'k','linewidth',2)
hold on
plot([0 1],[0 1],'k--')
plot(roc(disc_I,1),roc(disc_I,2),'r*','markersize',15,'linewidth',2);
text(roc(disc_I,1)+0.01,roc(disc_I,2)-0.05,'Sleep-wake cutoff','fontsize',15,'color','r');
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','southeast')
set(gca,'fontsize',15)

% B: Sleep surge histogram
all_pts_spikes_bins = sleep_hist_out.all_pts_spikes_bins;
all_pts_sleep_bins = sleep_hist_out.all_pts_sleep_bins;
sig_bins = sleep_hist_out.sig_bins;

nexttile
spike_bins = nanmean(all_pts_spikes_bins,1);
sleep_bins = nanmean(all_pts_sleep_bins,1);
times = linspace(-12,12,length(spike_bins));
spp = plot(times,spike_bins,'k-','linewidth',2);
hold on
plot(times(sig_bins),spike_bins(sig_bins)+0.05,'r*','linewidth',2);
hold on
%slp = plot(times,sleep_bins,'k--');
title('Spike rate surrounding sleep onset')
xlabel('Hours')
ylabel('Spikes/elecs/min')
xlim([-12 12])
yl = ylim;
plot([0 0],yl,'k--','linewidth',2)
set(gca,'fontsize',15)
%legend([spp;slp],{'Spikes','Proportion in sleep'})
ylim(yl);

% C: Sleep vs wake overall rate
nexttile
all_rates = bin_out.all_rates;
plot_paired_data(all_rates',{'Wake','Sleep'},'Spike/elec/min','paired','violin')
title('Spike rate in sleep vs wake')

% D: Sleep vs wake co-spiking
nexttile
all_coi = bin_out.all_coi;
plot_paired_data(all_coi',{'Wake','Sleep'},'Spike COI','paired','violin')
title('Co-spiking in sleep vs wake')

% E: Sleep vs wake node strength
nexttile
ns_sw = bin_out.ns_sw;
plot_paired_data(ns_sw',{'Wake','Sleep'},'Average node strength','paired','violin')
title('Functional connectivity in sleep vs wake')

print([out_folder,'Fig2'],'-dpng')
close(gcf)

%% Figure 3 - Effect of anatomical location
figure
set(gcf,'position',[10 10 1300 700])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

% A: Spike rate according to sleep vs wake and anatomical location
nexttile
main = bin_out.main;
r_ad_ana = bin_out.r_ad_ana;
interaction_plot_and_stats(r_ad_ana{1},main{1},'Spike/elec/min',{'Wake','Sleep'},0,'violin');
title('Spike rate by anatomical location')

% B: Spike rate according to sleep vs wake and SOZ
nexttile
r_ad_ana = bin_out.r_ad_ana;
interaction_plot_and_stats(r_ad_ana{3},main{3},'Spike/elec/min',{'Wake','Sleep'},0,'violin');
title('Spike rate by SOZ vs not SOZ')

% C: Spike latency according to sleep vs wake and anatomical location
nexttile
r_rl_ana = bin_out.r_rl_ana;
interaction_plot_and_stats(r_rl_ana{1}*1e3,main{1},'Spike latency (ms)',{'Wake','Sleep'},0,'violin');
title('Spike latency by anatomical location')

% D: Spike latency according to sleep vs wake and SOZ
nexttile
r_rl_ana = bin_out.r_rl_ana;
interaction_plot_and_stats(r_rl_ana{3}*1e3,main{3},'Spike latency (ms)',{'Wake','Sleep'},0,'violin');
title('Spike latency by SOZ vs not SOZ')

print([out_folder,'Fig3'],'-dpng')
close(gcf)

%% Figure 4 - Effect of seizures
figure
set(gcf,'position',[10 10 1400 700])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

sig_bins = sz_out.sig_bins;
all_pts_spikes_bins = sz_out.all_pts_spikes_bins;
all_pts_sleep_bins = sz_out.all_pts_sleep_bins;
surround_hours = sz_out.surround_hours;
sp_bins = nanmean(all_pts_spikes_bins,1)';
sleep_bins = nanmean(all_pts_sleep_bins,1)';
times = linspace(-surround_hours,surround_hours,length(sp_bins));

% A: Seizure surge histogram with superimposed sleep
nexttile
yLabels = {'Spikes/elec/min','Proportion asleep'};
h = stackedplot(times,[sp_bins,sleep_bins],'DisplayLabels',yLabels);
ax = findobj(h.NodeChildren, 'Type','Axes');
pause(0.1) % delete at your own risk
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
h.LineWidth = 2;
for k = 1:length(h.LineProperties)
    h.LineProperties(k).Color = 'k'; 
end
arrayfun(@(h)xline(h,0,'k--','LineWidth',2),ax) % plot line at 0
%plot(ax(1),times(sig_bins),sp_bins(sig_bins)+0.05,'r*','markersize',15,'linewidth',2); % sig times
set(gca,'fontsize',15)
title('Spikes and sleep surrounding seizures')

% B: Distribution amongst patients in post ictal change and pre-ictal
% change
nexttile
spikes_baseline = nanmean(all_pts_spikes_bins(:,baseline),2);
spikes_pre = nanmean(all_pts_spikes_bins(:,preictal),2);
spikes_post = nanmean(all_pts_spikes_bins(:,postictal),2);
plot_paired_data([spikes_baseline,spikes_pre,spikes_post]',...
    {'Baseline','Pre-ictal','Post-ictal'},'Spikes/elec/min','paired','violin')
title('Distribution in peri-ictal spike rates across patients')

% C: Post-ictal increase according to anatomical location
spikes_strat = sz_out.spikes_strat;
bl_post_strat = cell(length(spikes_strat),1);
for i = 1:length(bl_post_strat)
    npts = size(spikes_strat{i},2);
    ngroups = size(spikes_strat{i},1);
    nbins = size(spikes_strat{i},3);
    bl_post_strat{i} = nan(ngroups,npts,2); % bl and post
end
for g = 1:length(spikes_strat)
    bl_post_strat{g}(:,:,1) = nanmean(spikes_strat{g}(:,:,baseline),3);
    bl_post_strat{g}(:,:,2) = nanmean(spikes_strat{g}(:,:,postictal),3);
end

% anatomical localiztion
nexttile
interaction_plot_and_stats(bl_post_strat{1},main{1},'Spike/elec/min',...
    {'Baseline','Post-ictal'},0,'violin');
title('Post-ictal spike rate by anatomical location')

% D: Post-ictal increase according to SOZ vs not
nexttile
interaction_plot_and_stats(bl_post_strat{3},main{3},'Spike/elec/min',...
    {'Baseline','Post-ictal'},0,'violin');
title('Post-ictal spike rate by SOZ vs not SOZ')
print([out_folder,'Fig4'],'-dpng')
close(gcf)

end