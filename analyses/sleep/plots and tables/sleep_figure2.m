function sleep_figure2

%{
Note that the results text for the 3 and 12 hours peri-ictal windows are
hard coded in. IF I CHANGE THE ANALYSES I NEED TO CHANGE THE HARD CODING.
%}

%% Parameters
plot_type = 'scatter';
nblocks = 6;

%{
myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];
%}

myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];



locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
%out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');
fprintf(fid,'<p><b>Changes in spikes with seizures</b>');

figure
set(gcf,'position',[100 100 1200 1000])
tiledlayout(3,3,'tilespacing','tight','padding','tight')

%% Seizure time of day
nexttile([1 1])
all_tod_rate = sz_circ_out.all_tod_rate;
tod_edges = bin_out.tod_edges;
counts = nansum(all_tod_rate,1);
polar = convert_times_to_polar(tod_edges,'radians');

hours_mins = convert_edges_to_hours_mins(tod_edges);
nbins = length(counts);
skip = 18;

polarhistogram('BinEdges',polar,'BinCounts',counts,...
    'displayStyle','stairs','linewidth',2,'edgecolor',[0.2 0.2 0.2])
set(gca,'ThetaDir','clockwise');
set(gca,'ThetaZeroLocation','top');
thetaticks(polar(1:skip:nbins)*360/(2*pi))
thetaticklabels(hours_mins(1:skip:nbins+1))
set(gca,'fontsize',15)
title('Seizure count')
[pval z all_mu] = test_pt_circular_means(all_tod_rate,polar,hours_mins);

text(7*pi/8,15,get_p_text(pval),'fontsize',15,'horizontalalignment','center');
%{
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
%}

fprintf(fid,['<br>We asked whether patients demonstrated a circadian rhythm '...
    'in their seizure times. The distribution of seizure counts by time of day '...
    'was not significantly different from a uniform distribution (Rayleigh test '...
    'of patients'' circular means: z = %1.1f, %s) (Fig. 3A).'],...
    z,get_p_html(pval));


%% Percent asleep
nexttile([1 1])
sz_rate_sw = sz_circ_out.sz_rate_sw;
pause(0.3)
stats = plot_paired_data(sz_rate_sw',{'wake','sleep','sleep'},'Seizures/min','paired',plot_type);
pause(0.3)
xlim([0 0.11])
ylim([0 0.11])
title('Seizure frequency in wake and sleep')

% Results text
fprintf(fid,[' Seizure rates were not significantly different between sleep (median %1.3f seizures/min)'...
    ' and wake (median %1.3f seizures/min) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 3B).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));



%{

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
%}

nexttile([1 1])
pre_wake = sz_circ_out.pre_wake;
n_sleep_wake = bin_out.n_sleep_wake;
perc_sz_asleep = cellfun(@(x) prc_asleep(x),pre_wake);

loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
[p,~,stats] = ranksum(perc_sz_asleep(temporal),perc_sz_asleep(extra));
plot(1+randn(sum(temporal),1)*0.05,perc_sz_asleep(temporal),'o','linewidth',2,'color',myColours(2,:))
hold on
plot([0.7 1.3],[nanmedian(perc_sz_asleep(temporal)) nanmedian(perc_sz_asleep(temporal))],...
    'linewidth',2,'color',myColours(2,:))
plot(2+randn(sum(extra),1)*0.05,perc_sz_asleep(extra),'o','linewidth',2,'color',myColours(3,:))
plot([1.7 2.3],[nanmedian(perc_sz_asleep(extra)) nanmedian(perc_sz_asleep(extra))],...
    'linewidth',2,'color',myColours(3,:))
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})

ylabel('Seizures arising from sleep (%)')
title('Seizures arising from sleep')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)
yt = yticks;
ytl = get(gca,'yticklabels');
%yticks(yt(1:end-1));
yticklabels(ytl(1:end-1))

W = stats.ranksum;
nt = sum(~(isnan(perc_sz_asleep(temporal))));
ne = sum(~(isnan(perc_sz_asleep(extra))));
%ns = min([sum(temporal),sum(extra)]);
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,[' The proportion of seizures arising from sleep was similar '...
    'between patients with temporal lobe epilepsy (median = %1.1f%%) and extra-temporal'...
    ' lobe epilepsy (median = %1.1f%%) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>temporal</sub></i> = %d, <i>N<sub>extra-temporal</sub></i> = %d) ='...
    ' %1.1f, %s) (Fig. 3C).'],nanmedian(perc_sz_asleep(temporal)),nanmedian(perc_sz_asleep(extra)),...
    nt,ne,U,get_p_html(p));

%% Supplemental
%{
mt = contains(loc,'mesial temporal');
nc = contains(loc,'cortical') | contains(loc,'other cortex'); % temporal neocortical and other cortex

[p,~,stats] = ranksum(perc_sz_asleep(mt),perc_sz_asleep(nc));
W = stats.ranksum;
nt = sum(~(isnan(perc_sz_asleep(mt))));
ne = sum(~(isnan(perc_sz_asleep(nc))));
%ns = min([sum(temporal),sum(extra)]);
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,[' The proportion of seizures arising from sleep was also similar '...
    'between patients with mesial temporal epilepsy (median = %1.1f%%) and neocortical'...
    ' epilepsy (median = %1.1f%%) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>mesial temporal</sub></i> = %d, <i>N<sub>neocortical</sub></i> = %d) ='...
    ' %1.1f, %s) (Supplementary Fig. 1C).</p>'],nanmedian(perc_sz_asleep(mt)),nanmedian(perc_sz_asleep(nc)),...
    nt,ne,U,get_p_html(p));
%}

fprintf(fid,[' The proportion of seizures arising from sleep was also similar '...
    'between patients with mesial temporal epilepsy and neocortical epilepsy '...
    '(Supplementary Fig. 1C, Supplementary Table 2).</p>']);


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
shaded_error_bars_fc(times,median_spikes,iqr_spikes,myColours(1,:));
hold on
plot([0 0],ylim,'k--','linewidth',3)
ylabel('Spikes/elecs/min')
title('Peri-ictal spike rates and % asleep')
set(gca,'fontsize',15)
%xlabel('Hours surrounding seizure')

fprintf(fid,['<p>We next looked at spike rates surrounding seizures. '...
    'Visually, there was an increase in spike rates after seizures that '...
    'tracked with a postictal increase in sleep classification (Fig. 3D).']);

%% Pre vs post ictal
nexttile([1 1])
nbins = size(all_pts_spikes_bins,2);
pre = 1:nbins/2;
post = nbins/2+1:nbins;
pre_post = [nanmean(all_pts_spikes_bins(:,pre),2),nanmean(all_pts_spikes_bins(:,post),2)];
stats = plot_paired_data(pre_post',{'preictal','postictal','postictally'},'Spikes/elec/min','paired',plot_type);
title('Pre- vs postictal spike rates')

% Results text
fprintf(fid,[' Across all patients, the median spike rate postictally (median %1.2f spikes/elecs/min)'...
    ' was higher than that preictally (median %1.2f spikes/elecs/min) '...
    '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 3E).'...
    ' This was also true when the pre- and postictal windows were defined to be 3 hours each '...
    'and when they were defined to be 12 hours each (p < 0.001 for each of the '...
    'alternative peri-ictal time windows) (Supplemental Table 1).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));


%% Temporal vs extratemporal
nexttile([1 1])
loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
rate_diff = (pre_post(:,2)-pre_post(:,1));
[p,~,stats] = ranksum(rate_diff(temporal),rate_diff(extra));
plot(1+randn(sum(temporal),1)*0.05,rate_diff(temporal),'o','linewidth',2,'color',myColours(2,:))
hold on
plot([0.7 1.3],[nanmedian(rate_diff(temporal)) nanmedian(rate_diff(temporal))],...
    'linewidth',2,'color',myColours(2,:))
plot(2+randn(sum(extra),1)*0.05,rate_diff(extra),'o','linewidth',2,'color',myColours(3,:))
plot([1.7 2.3],[nanmedian(rate_diff(extra)) nanmedian(rate_diff(extra))],...
    'linewidth',2,'color',myColours(3,:))
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})
ylabel('Post-pre spikes/elec/min')
title('Post-preictal spike rate difference')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)

% text
W = stats.ranksum;
nt = sum(~(isnan(rate_diff(temporal))));
ne = sum(~(isnan(rate_diff(extra))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,[' Patients with temporal lobe epilepsy had a greater '...
    'pre-to-postictal increase in spike rates '...
    '(median = %1.2f spikes/elecs/min) relative to those with extra-temporal'...
    ' lobe epilepsy (median = %1.2f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>temporal</sub></i> = %d, <i>N<sub>extra-temporal</sub></i> = %d) ='...
    ' %1.1f, %s) (Fig. 3F). This was also true when the pre- and post-ictal'...
    ' windows were defined to be 3 hours each (<i>p</i> = 0.035) and when they '...
    'were defined to be 12 hours each (<i>p</i> < 0.001) (Supplemental Table 1).'],nanmedian(rate_diff(temporal)),nanmedian(rate_diff(extra)),...
    nt,ne,U,get_p_html(p));

%% Supplemental
%{
[p,~,stats] = ranksum(rate_diff(mt),rate_diff(nc));
% text
W = stats.ranksum;
nt = sum(~(isnan(rate_diff(mt))));
ne = sum(~(isnan(rate_diff(nc))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,[' Patients with mesial temporal epilepsy also had a greater '...
    'pre-to-postictal increase in spike rates '...
    '(median = %1.2f spikes/elecs/min) than did those with neocortical'...
    ' epilepsy (median = %1.2f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>mesial temporal</sub></i> = %d, <i>N<sub>neocortical</sub></i> = %d) ='...
    ' %1.1f, %s) (Supplementary Fig. 1D).</p>'],nanmedian(rate_diff(mt)),nanmedian(rate_diff(nc)),...
    nt,ne,U,get_p_html(p));
%}
fprintf(fid,[' Patients with mesial temporal epilepsy also had a greater '...
    'pre-to-postictal increase in spike rates '...
    'than did those with neocortical'...
    ' epilepsy (Supplementary Fig. 1D, Supplementary Table 2).</p>']);

%% Percent detected asleep
nexttile([1 1])
plot(times,mean_sleep,'color',[0.9290, 0.6940, 0.1250],'linewidth',2)
hold on

set(gca,'fontsize',15)
ylabel('% Classified asleep')
%title('Peri-ictal sleep classification')
xlabel('Hours surrounding seizure')
plot([0 0],ylim,'k--','linewidth',3)

%% Pre-ictal
nbins = size(all_pts_spikes_bins,2);
early_pre = 1:nbins/4; % first quarter
late_pre = nbins/4+1:nbins/2; % second quarter
early_late = [nanmean(all_pts_spikes_bins(:,early_pre),2),nanmean(all_pts_spikes_bins(:,late_pre),2)];

nexttile
% Early vs late preictal spike rates
stats = plot_paired_data(early_late',{'early preictal','late preictal','late'},'Spikes/elec/min','paired',plot_type);
title('Early vs late preictal spike rates')

% Results text
fprintf(fid,['<p>Visually, there was no obvious preictal increase in spikes (Fig. 3D).'...
    ' Across all patients, the median spike rate in the late preictal period (median %1.2f spikes/elecs/min)'...
    ' was similar to that of the early preictal period (median %1.2f spikes/elecs/min) when examining 6-hour preictal periods '...
    '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 3G).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

% Secondary analyses
fprintf(fid,[' However, spike rates were significantly higher in the late preictal period than'...
    ' in the early preictal period when performing our secondary analysis of 3-hour preictal windows (<i>p</i> = 0.018) '...
    'and 12-hour preictal windows (<i>p</i> = 0.020) (Supplemental Table 1).']); 

%% Rate diff in preictal period by epilepsy localization
nexttile
rate_diff = early_late(:,2)-early_late(:,1);
[p,~,stats] = ranksum(rate_diff(temporal),rate_diff(extra));
plot(1+randn(sum(temporal),1)*0.05,rate_diff(temporal),'o','linewidth',2,'color',myColours(2,:))
hold on
plot([0.7 1.3],[nanmedian(rate_diff(temporal)) nanmedian(rate_diff(temporal))],...
    'linewidth',2,'color',myColours(2,:))
plot(2+randn(sum(extra),1)*0.05,rate_diff(extra),'o','linewidth',2,'color',myColours(3,:))
plot([1.7 2.3],[nanmedian(rate_diff(extra)) nanmedian(rate_diff(extra))],...
    'linewidth',2,'color',myColours(3,:))
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})
ylabel('Late-early spikes/elec/min')
title('Late-early preictal spike rate difference')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)

% text
W = stats.ranksum;
nt = sum(~(isnan(rate_diff(temporal))));
ne = sum(~(isnan(rate_diff(extra))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);

fprintf(fid,[' For the primary analysis of a 6-hour preictal period, the early-to-late'...
    ' preictal spike rate change was similar between patients with temporal lobe epilepsy '...
    '(median = %1.2f spikes/elecs/min) and those with extra-temporal lobe epilepsy'...
    ' (median = %1.2f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>temporal</sub></i> = %d, <i>N<sub>extra-temporal</sub></i> = %d) ='...
    ' %1.1f, %s) (Fig. 3H). There was also no significant difference in the preictal spike '...
    'rate change between epilepsy localizations when studying the 3-hour '...
    '(<i>p</i> = 0.66) or 12-hour preictal windows (<i>p</i> = 0.42) (Supplemental Table 1).'],nanmedian(rate_diff(temporal)),nanmedian(rate_diff(extra)),...
    nt,ne,U,get_p_html(p));

%% Supplemental
%{
[p,~,stats] = ranksum(rate_diff(mt),rate_diff(nc));
W = stats.ranksum;
nt = sum(~(isnan(rate_diff(mt))));
ne = sum(~(isnan(rate_diff(nc))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);

fprintf(fid,[' The early-to-late'...
    ' preictal spike rate change was also similar between patients with mesial temporal epilepsy '...
    '(median = %1.2f spikes/elecs/min) and those with neocortical epilepsy'...
    ' (median = %1.2f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>mesial temporal</sub></i> = %d, <i>N<sub>neocortical</sub></i> = %d) ='...
    ' %1.1f, %s) (Supplementary Fig. 1E). Overall, these results imply that there is no '...
    'clear preictal spike rate change, and the preictal spike rate change does'...
    ' not vary by epilepsy localization. </p>'],nanmedian(rate_diff(mt)),nanmedian(rate_diff(nc)),...
    nt,ne,U,get_p_html(p));
%}
fprintf(fid,[' The early-to-late'...
    ' preictal spike rate change was also similar between patients with mesial temporal epilepsy '...
    'and those with neocortical epilepsy (Supplementary Fig. 1E, Supplementary Table 2).'...
    ' Overall, these results imply that there is no '...
    'consistent preictal spike rate change, and the preictal spike rate change does'...
    ' not vary by epilepsy localization. </p>']);
%% Add annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.33 0.9 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0.67 0.9 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.545 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0.33 0.545 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0.67 0.545 0.1 0.1],'String','F','fontsize',25,'linestyle','none')
%annotation('textbox',[0 0.23 0.1 0.1],'String','G','fontsize',25,'linestyle','none')
annotation('textbox',[0.33 0.23 0.1 0.1],'String','G','fontsize',25,'linestyle','none')
annotation('textbox',[0.67 0.23 0.1 0.1],'String','H','fontsize',25,'linestyle','none')

%% Add box around D
annotation('rectangle',[0.001 0.001 0.32 0.65],'color','k','linewidth',2,'linestyle','--')

print([out_folder,'Fig3'],'-depsc')
close(gcf)

%% Bonus analysis (probably supplemental figure) looking for pre-ictal spike change
%{
figure
set(gcf,'position',[100 100 1200 334]);
tiledlayout(1,3,'tilespacing','tight','padding','tight')


nexttile
% Early vs late preictal sleep
early_late_sleep = [nanmean(all_pts_sleep_bins(:,early_pre),2),nanmean(all_pts_sleep_bins(:,late_pre),2)];
stats = plot_paired_data(early_late_sleep',{'early preictal period','late preictal period','late preictal period'},'Proportion asleep','paired',plot_type);
title('Early vs late preictal sleep')

    
    
    
    There was a concomitant increase '...
    'in the proportion of patients detected to be asleep preictally when studying the '...
    '3-hour preictal window ('...
    'early median 0.22 spikes/elecs/min, late median 0.23 spikes/elecs/min,'...
    ' <i>T<sup>+</sup></i> = 1199, <i>p</i> = 0.043) but not when studying the 12-hour preictal window '...
    '(early median 0.23 spikes/elecs/min, late median 0.23 spikes/elecs/min,'...
    ' <i>T<sup>+</sup></i> = 1492.5, <i>p</i> = 0.074). This suggests the possibility that some of the '...
    'apparent preictal increase in spike rates is related to a preictal tendency to '...
    'be asleep.CHECK ME']);
    %}



%print([out_folder,'SuppFig1'],'-depsc')



end

function prc = prc_asleep(x)

prc = 100 * sum(x==0)/(sum(x==1)+sum(x==0));

end