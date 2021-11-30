function sleep_plots(out,do_save)

%% Parameters
plot_type = 'scatter';
nblocks = 6;
myColours = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250];

%% Unpack substructures
unpack_any_struct(out);

%% Fig 1 - circadian analysis
figure
set(gcf,'position',[100 100 1400 600])
tiledlayout(2,3,'tilespacing','tight','padding','tight')

%% A - PSD
median_psd = circ_out.median_psd;
iqr_psd = circ_out.iqr_psd;
periods = circ_out.periods;

nexttile
shaded_error_bars(periods,median_psd,iqr_psd,[]);
xlim([0 100])
xlabel('Period (hours)')
ylabel({'Spike rate', 'normalized power spectrum'});
set(gca,'fontsize',15)
title('Spike rate power spectral density')

%% B: ROC curve
roc = roc_out.roc;
auc = roc_out.auc;
disc_I = roc_out.disc_I;
nexttile
plot(roc(:,1),roc(:,2),'linewidth',2)
hold on
plot([0 1],[0 1],'k--')
plot(roc(disc_I,1),roc(disc_I,2),'*','markersize',15,'linewidth',2,'color',myColours(2,:));
text(roc(disc_I,1)+0.01,roc(disc_I,2)-0.05,'Sleep-wake cutoff','fontsize',15,'color',myColours(2,:));
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','southeast')
set(gca,'fontsize',15)

%% C: Sleep surge histogram
all_pts_spikes_bins = sleep_hist_out.all_pts_spikes_bins;
all_pts_sleep_bins = sleep_hist_out.all_pts_sleep_bins;
sig_bins = sleep_hist_out.sig_bins;

nexttile
median_spikes = nanmedian(all_pts_spikes_bins,1);
mean_sleep = nanmean(all_pts_sleep_bins,1)*100;
iqr_spikes = prctile(all_pts_spikes_bins,[25,75],1);
times = linspace(-12,12,length(median_spikes));
ylabels = ["Spike rate","% asleep"];
s = stackedplot(times,[median_spikes',mean_sleep'],'linewidth',2,...
    "DisplayLabels",ylabels);

for k = 1:length(s.LineProperties)
    if k <= size(myColours,1)
        s.LineProperties(k).Color = myColours(k,:);
    end
end
ax = findobj(s.NodeChildren, 'Type','Axes');
arrayfun(@(h)xline(h,0,'--k','LineWidth',2),ax)
pause(0.3)
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
pause(0.3)
xlabel('Hours')
xlim([-12 12])
set(gca,'fontsize',15)
%title('Spike rate surrounding sleep onset')


%% D-F: sleep pca stuff
sleep_bins = out.sleep_hist_out.all_pts_spikes_bins;
% Subtract mean
sleep_bins = (sleep_bins - nanmean(sleep_bins,2))./nanstd(sleep_bins,[],2);
[coeff,score,latent] = pca(sleep_bins,'Rows','pairwise');
locs = out.circ_out.all_locs;
% Top scorers
[~,top] = max(score,[],1);
% Bottom scorers
[~,bottom] = min(score,[],1);
top_12 = [top(1:2),bottom(1:2)];

% D: Variances
nexttile
stem(latent,'linewidth',2)
title('Principal component variances')
set(gca,'fontsize',15)
xlabel('Component')
ylabel('Variance')

% E: First two components
nexttile
p1= plot(times,coeff(:,1),'linewidth',2);
hold on
p2 = plot(times,coeff(:,2),'linewidth',2);
plot([0 0],ylim,'k--','linewidth',2)
%p3 = plot(times,coeff(:,3),'linewidth',2);
ylabel('Coefficients')
title('Top principal components')
legend([p1 p2],{'Component 1','Component 2'},'fontsize',15,'location','northwest')
set(gca,'fontsize',15)
xlabel('Hours')
xlim([times(1) times(end)])

% F: Localization
nexttile
c = 1;
tloc = strcmp(locs,'temporal');
oloc = strcmp(locs,'other');
plot(1+0.05*randn(sum(tloc),1),score(tloc,c),'o','linewidth',2)
hold on
plot(2+0.05*randn(sum(oloc),1),score(oloc,c),'o','linewidth',2)
xlim([0.5 2.5])
plot(xlim,[0 0],'k--','linewidth',2)
p = ranksum(score(tloc,c),score(oloc,c));
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
newyl = [yl(1) yl(1) + 1.17*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'horizontalalignment','center','fontsize',15)
title(sprintf('Component %d score by localization',c))
set(gca,'fontsize',15)
ylim(newyl)
xticks([1 2])
ylabel('Score')
xticklabels({'Temporal','Extra-temporal'})

if do_save
print([out_folder,'Fig1'],'-dpng')
close(gcf)
end

%% Fig 2 - exploring sleep analysis
figure
set(gcf,'position',[100 100 1400 600])
tiledlayout(2,3,'tilespacing','tight','padding','tight')

%% 2A overall spike rates
all_rate = bin_out.all_rates;
nexttile
plot_paired_data(all_rate(:,1:2)',{'wake','sleep'},'Spikes/elec/min','paired',plot_type)
title('Overall spike rates')

%% 2B independent spikes
seq_sw = bin_out.seq_sw;
nexttile
plot_paired_data(seq_sw(:,1:2)',{'wake','sleep'},'Spikes/elec/min','paired',plot_type)
title('# Independent spikes sources')

%% 2C spike spread
nexttile
plot_paired_data(seq_sw(:,3:4)',{'wake','sleep'},'# Spikes/sequence','paired',plot_type)
title('Spike spread')

%% 2D RL sleep-wake correlation across patients
rl_sw_corr = bin_out.rl_sw_corr;
nexttile
plot(rl_sw_corr,'o','linewidth',2)
hold on
plot(xlim,[0 0],'k--','linewidth',2)
ylabel('Correlation coefficient')
xlabel('Patient')
xticklabels([])
ylim([-1 1])
title('Sleep vs wake correlation in spike spread')
set(gca,'fontsize',15)

%% 2E NS
ns_sw = bin_out.ns_sw;
nexttile
plot_paired_data(ns_sw',{'wake','sleep'},'Average node strength','paired',plot_type)

%% 2F SOZ spike rate ranking
nexttile
soz_rank_sw = bin_out.soz_rank_sw;
plot_paired_data(soz_rank_sw',{'wake','sleep'},'Rank in spike rate','paired',plot_type,'ranking')
title({'SOZ spike rate ranking','by wake vs sleep'})

if do_save
print([out_folder,'Fig2'],'-dpng')
close(gcf)
end

%% Figure 3 - spikes tend to go up after seizures, depends on localization
figure
set(gcf,'position',[100 100 900 700])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

%sig_bins = sz_out.sig_bins;
all_pts_spikes_bins = sz_out.all_pts_spikes_bins;
all_pts_sleep_bins = sz_out.all_pts_sleep_bins;
surround_hours = sz_out.surround_hours;
sp_bins = nanmean(all_pts_spikes_bins,1)';
sleep_bins = nanmean(all_pts_sleep_bins,1)';
times = linspace(-surround_hours,surround_hours,length(sp_bins));

%% 3A: Seizure surge histogram with superimposed sleep
nexttile
yLabels = {'Spikes/elec/min','Proportion asleep'};
h = stackedplot(times,[sp_bins,sleep_bins],'DisplayLabels',yLabels);
ax = findobj(h.NodeChildren, 'Type','Axes');
pause(0.1) % delete at your own risk
set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
h.LineWidth = 2;
for k = 1:length(h.LineProperties)
    if k <= size(myColours,1)
        h.LineProperties(k).Color = myColours(k,:);
    end
end
arrayfun(@(h)xline(h,0,'k--','LineWidth',2),ax) % plot line at 0
%plot(ax(1),times(sig_bins),sp_bins(sig_bins)+0.05,'r*','markersize',15,'linewidth',2); % sig times
set(gca,'fontsize',15)
xlabel('Hours')
title('Spikes and sleep surrounding seizures')

%% PCA stuff
% Subtract mean
sp_bins = all_pts_spikes_bins;
sp_bins = (sp_bins - nanmean(sp_bins,2))./nanstd(sp_bins,[],2);
[coeff,score,latent] = pca(sp_bins,'Rows','complete');
locs = out.circ_out.all_locs;
% Top scorers
[~,top] = max(score,[],1);
% Bottom scorers
[~,bottom] = min(score,[],1);
top_12 = [top(1:2),bottom(1:2)];

%% 3B: Variances
nexttile
stem(latent,'linewidth',2)
title('Principal component variances')
set(gca,'fontsize',15)
xlabel('Component')
ylabel('Variance')

%% 3C: First two principal components
nexttile
p1= plot(times,coeff(:,1),'linewidth',2);
hold on
p2 = plot(times,coeff(:,2),'linewidth',2);
plot([0 0],ylim,'k--','linewidth',2)
%p3 = plot(times,coeff(:,3),'linewidth',2);
ylabel('Coefficients')
title('Top principal components')
legend([p1 p2],{'Component 1','Component 2'},'fontsize',15,'location','northwest')
set(gca,'fontsize',15)
xlabel('Hours')
xlim([times(1) times(end)])

%% 3D: First PCA score based on sz localization
nexttile
c = 1;
tloc = strcmp(locs,'temporal');
oloc = strcmp(locs,'other');
plot(1+0.05*randn(sum(tloc),1),score(tloc,c),'o','linewidth',2)
hold on
plot(2+0.05*randn(sum(oloc),1),score(oloc,c),'o','linewidth',2)
xlim([0.5 2.5])
plot(xlim,[0 0],'k--','linewidth',2)
p = ranksum(score(tloc,c),score(oloc,c));
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
newyl = [yl(1) yl(1) + 1.17*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'horizontalalignment','center','fontsize',15)
title(sprintf('Component %d score by localization',c))
set(gca,'fontsize',15)
ylim(newyl)
xticks([1 2])
ylabel('Score')
xticklabels({'Temporal','Extra-temporal'})

if do_save
print([out_folder,'Fig3'],'-dpng')
close(gcf)
end

%{
%% Bonus test - SOZ ranking sw
figure
soz_rank_sw = bin_out.soz_rank_sw;
plot_paired_data(soz_rank_sw',{'Wake','Sleep'},'Rank in spike rate','paired',plot_type)
title({'SOZ rank in spike rate','by wake vs sleep'})
print([out_folder,'SOZTest'],'-dpng')
close(gcf)

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
plot_paired_data(circ_P{1},main_locs,'Relative circadian power','paired',plot_type)
title({'Relative circadian power','by spike location'})

% Compare cyclical power between patients with different soz localizations
nexttile
mt = cellfun(@(x) strcmp(x,'temporal'), all_locs);
circ_mt = all_circ_P(mt);
circ_other = all_circ_P(~mt);
circ_mt = [circ_mt;nan(length(mt)-length(circ_mt),1)]; % pad both with nans
circ_other = [circ_other;nan(length(mt)-length(circ_other),1)];
plot_paired_data(([circ_mt,circ_other])',{'Temporal','Other'},'Relative circadian power',...
    'unpaired',plot_type)
title({'Relative circadian power','by SOZ localization'});
print([out_folder,'Fig1'],'-dpng')
close(gcf)

%% Figure 2 - What happens to spikes with sleep
figure
set(gcf,'position',[10 10 1000 1000])
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
median_spikes = nanmedian(all_pts_spikes_bins,1);
iqr_spikes = prctile(all_pts_spikes_bins,[25,75],1);
times = linspace(-12,12,length(median_spikes));
%spp = plot(times,spike_bins,'k-','linewidth',2);
shaded_error_bars(times,median_spikes,iqr_spikes,[0 0 0]);
hold on
plot(times(sig_bins),median_spikes(sig_bins)+0.05,'r*','linewidth',2);
hold on
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
%{
nexttile
all_rates = bin_out.all_rates;
plot_paired_data(all_rates',{'Wake','Sleep'},'Spike/elec/min','paired',plot_type)
title('Spike rate in sleep vs wake')
%}


nseq = bin_out.seq_sw(:,1:2);
seq_length = bin_out.seq_sw(:,3:4);

% C: sleep vs wake nseq
nexttile
plot_paired_data(nseq',{'Wake','Sleep'},'Number of sequences','paired',plot_type)

% D sleep vs wake sequence length
nexttile
plot_paired_data(seq_length',{'Wake','Sleep'},'Number of spikes in sequences','paired',plot_type)

% D: Sleep vs wake co-spiking
%{
nexttile
all_coi = bin_out.all_coi;
plot_paired_data(all_coi',{'Wake','Sleep'},'Spike COI','paired',plot_type)
title('Co-spiking in sleep vs wake')
%}

% E: Sleep vs wake node strength
nexttile
ns_sw = bin_out.ns_sw;
plot_paired_data(ns_sw',{'Wake','Sleep'},'Average node strength','paired',plot_type)
title('Functional connectivity in sleep vs wake')

print([out_folder,'Fig2'],'-dpng')
close(gcf)

%% Figure 3 - Effect of anatomical location
figure
set(gcf,'position',[10 10 1400 900])
tiledlayout(2,4,'tilespacing','tight','padding','tight')
main = bin_out.main;
r_ad_ana = bin_out.r_ad_ana;

%
% A: Spike rate according to sleep vs wake and anatomical location
nexttile
interaction_plot_and_stats(r_ad_ana{1},make_multi_line(main{1}),'Spike/elec/min',{'Wake','Sleep'},0,plot_type);
title('Spike rate by anatomical location')

% Spike rate according to sleep vs wake and SOZ
nexttile
r_ad_ana = bin_out.r_ad_ana;
interaction_plot_and_stats(r_ad_ana{3},main{3},'Spike/elec/min',{'Wake','Sleep'},0,plot_type);
title('Spike rate by SOZ vs not SOZ')

% C: Spike latency according to sleep vs wake and anatomical location
nexttile
r_rl_ana = bin_out.r_rl_ana;
interaction_plot_and_stats(r_rl_ana{1}*1e3,make_multi_line(main{1}),'Spike latency (ms)',{'Wake','Sleep'},0,plot_type);
title('Spike latency by anatomical location')

% Spike latency according to sleep vs wake and SOZ
nexttile
r_rl_ana = bin_out.r_rl_ana;
interaction_plot_and_stats(r_rl_ana{3}*1e3,main{3},'Spike latency (ms)',{'Wake','Sleep'},0,plot_type);
title('Spike latency by SOZ vs not SOZ')


% Is sleep-related increase in spike rate difference across anatomical
% locations
nexttile
plot_and_stats_change(r_ad_ana{1},main{1},'Spike/elec/min','paired')
title('Sleep-related rate change')

% Is sleep-related increase in spike rate higher for SOZ?
nexttile
plot_and_stats_change(r_ad_ana{3},main{3},'Spike/elec/min','paired')
title('Sleep-related rate change')


% Spike latency increase based on anatomical location
nexttile
plot_and_stats_change(r_rl_ana{1},main{1},{'Latency difference (ms)'},'paired')
title('Sleep-related latency change')


% Is sleep-related increase in latency higher for SOZ?
nexttile
plot_and_stats_change(r_rl_ana{3},main{3},{'Latency difference (ms)'},'paired')
title('Sleep-related latency change')

print([out_folder,'Fig3'],'-dpng')
close(gcf)

%% Figure 4 - Effect of seizures
figure
set(gcf,'position',[10 10 1300 800])
tiledlayout(2,3,'tilespacing','tight','padding','tight')

sig_bins = sz_out.sig_bins;
all_pts_spikes_bins = sz_out.all_pts_spikes_bins;
all_pts_sleep_bins = sz_out.all_pts_sleep_bins;
surround_hours = sz_out.surround_hours;
sp_bins = nanmean(all_pts_spikes_bins,1)';
sleep_bins = nanmean(all_pts_sleep_bins,1)';
times = linspace(-surround_hours,surround_hours,length(sp_bins));

% Seizure surge histogram
nexttile
median_spikes = nanmedian(all_pts_spikes_bins,1);
iqr_spikes = prctile(all_pts_spikes_bins,[25 75],1);
shaded_error_bars(times,median_spikes,iqr_spikes,[0 0 0]);
hold on
plot(times(sig_bins),sp_bins(sig_bins)+0.05,'r*','markersize',15,'linewidth',2); % sig times
plot([0 0],ylim,'k--','LineWidth',2);
set(gca,'fontsize',15)
title('Spikes surrounding seizures')
xlabel('Hours')
ylabel('Spikes/elec/min')

% A: Seizure surge histogram with superimposed sleep
%{
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
%}

% B: Distribution amongst patients in post ictal change and pre-ictal
% change
%{
nexttile
spikes_baseline = nanmean(all_pts_spikes_bins(:,baseline),2);
spikes_pre = nanmean(all_pts_spikes_bins(:,preictal),2);
spikes_post = nanmean(all_pts_spikes_bins(:,postictal),2);
plot_paired_data([spikes_baseline,spikes_pre,spikes_post]',...
    {'Baseline','Pre-ictal','Post-ictal'},'Spikes/elec/min','paired',plot_type)
title('Distribution in peri-ictal spike rates across patients')
%}

% C: Post-ictal increase according to anatomical location
spikes_strat = sz_out.spikes_strat;
bl_post_strat = cell(length(spikes_strat),1);
midpoint = size(spikes_strat{1},3)/2;
baseline = 1:nblocks; % one hour
preictal = midpoint-nblocks:midpoint;
postictal = midpoint+1:midpoint+nblocks;
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
interaction_plot_and_stats(bl_post_strat{1},make_multi_line(main{1}),'Spike/elec/min',...
    {'Baseline','Post-ictal'},0,plot_type);
title('Post-ictal spike rate by anatomical location')

% Post-ictal increase according to SOZ vs not
nexttile
interaction_plot_and_stats(bl_post_strat{3},main{3},'Spike/elec/min',...
    {'Baseline','Post-ictal'},0,plot_type);
title('Post-ictal spike rate by SOZ vs not SOZ')

% Seizure surge histogram for sleep
nexttile
median_sleep = nanmedian(all_pts_sleep_bins,1);
iqr_sleep = prctile(all_pts_sleep_bins,[25 75],1);
shaded_error_bars(times,median_sleep,iqr_sleep,[0 0 0]);
hold on
plot([0 0],ylim,'k--','LineWidth',2);
set(gca,'fontsize',15)
title('Sleep surrounding seizures')
xlabel('Hours')
ylabel('Proportion asleep')

% anatomical localiztion
nexttile
plot_and_stats_change(bl_post_strat{1},main{1},'Spike/elec/min',...
    'paired');
title('Post-ictal rate change')

% Post-ictal increase according to SOZ vs not
nexttile
plot_and_stats_change(bl_post_strat{3},main{3},'Spike/elec/min',...
    'paired');
title('Post-ictal rate change')


print([out_folder,'Fig4'],'-dpng')
close(gcf)
%}
end