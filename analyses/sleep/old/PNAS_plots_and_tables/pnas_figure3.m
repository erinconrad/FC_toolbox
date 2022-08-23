function pnas_figure3

%{
Note that the results text for the 3 and 12 hours peri-ictal windows are
hard coded in. IF I CHANGE THE ANALYSES I NEED TO CHANGE THE HARD CODING.
%}

%% Parameters
plot_type = 'scatter';

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
out_folder = [results_folder,'analysis/sleep/pnas/'];

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');
sfid = fopen([out_folder,'supplemental_results.html'],'a');
fprintf(fid,'<p><b>Changes in spikes with seizures</b>');

figure
set(gcf,'position',[100 100 1200 350])
tiledlayout(2,3,'tilespacing','tight','padding','tight')


%% Seizure time of day
all_tod_rate = sz_circ_out.all_tod_rate;
tod_edges = bin_out.tod_edges;
counts = nansum(all_tod_rate,1);
polar = convert_times_to_polar(tod_edges,'radians');

hours_mins = convert_edges_to_hours_mins(tod_edges);
nbins = length(counts);
[pval,z,all_mu] = test_pt_circular_means(all_tod_rate,polar,hours_mins);

fprintf(fid,['<br>We asked whether patients demonstrated a circadian rhythm '...
    'in seizure times. The distribution of seizure counts by time of day '...
    'was significantly different from a uniform distribution (Rayleigh test '...
    'of patients'' circular means: z = %1.1f, %s), with visual analysis '...
    'demonstrating more seizures in the morning and afternoon (Supplemental Figure 3A).'],...
    z,get_p_html(pval));

%% Percent asleep
sz_rate_sw = sz_circ_out.sz_rate_sw;
sz_rate_sw = sz_rate_sw/10; % seizures per minute (was previously seizures per ten minutes)
sz_rate_sw = sz_rate_sw*60*24; % seizures per day
stats = plot_paired_data(sz_rate_sw',{'wake','sleep','sleep'},'Seizures/day','paired',plot_type,0,0);

fprintf(fid,[' Seizure rates were similar between sleep (median %1.3f seizures/day)'...
    ' and wake (median %1.3f seizures/day) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Supplemental Figure 3B).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));


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
shaded_error_bars_fc(times,median_spikes,iqr_spikes,myColours(1,:));
hold on
plot([0 0],ylim,'k--','linewidth',3)
ylabel('Spikes/elecs/min')
title('Peri-ictal spike rates and % asleep')
set(gca,'fontsize',15)
%xlabel('Hours surrounding seizure')

fprintf(fid,['<p>We next examined spike rates surrounding seizures. '...
    'Visually, there was an increase in spike rates after seizures that '...
    'tracked with a postictal increase in sleep classification (Fig. 3A).']);

%% Pre vs post ictal
nexttile([2 1])
nbins = size(all_pts_spikes_bins,2);
pre = 1:nbins/2;
post = nbins/2+1:nbins;
pre_post = [nanmean(all_pts_spikes_bins(:,pre),2),nanmean(all_pts_spikes_bins(:,post),2)];
stats = plot_paired_data(pre_post',{'preictal','postictal','postictally'},'Spikes/elec/min','paired',plot_type);
title('Pre- vs postictal spike rates')

% Results text
fprintf(fid,[' Across all patients, the median spike rate postictally (median %1.2f spikes/elecs/min)'...
    ' was higher than that preictally (median %1.2f spikes/elecs/min) '...
    '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 3B).'...
    ' This was also true when the pre- and postictal windows were defined to be 3 hours each '...
    'and when they were defined to be 12 hours each (p < 0.001 for each of the '...
    'alternative peri-ictal time windows) (Supplemental Table 2).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));




%% Temporal vs extratemporal
nexttile([2 1])
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
    ' %1.1f, %s) (Fig. 3C). This was also true when the pre- and post-ictal'...
    ' windows were defined to be 3 hours each (<i>p</i> = 0.023) and when they '...
    'were defined to be 12 hours each (<i>p</i> < 0.001) (Supplemental Table 2).'],nanmedian(rate_diff(temporal)),nanmedian(rate_diff(extra)),...
    nt,ne,U,get_p_html(p));

fprintf(fid,[' Patients with mesial temporal epilepsy also had a greater '...
    'pre-to-postictal increase in spike rates '...
    'than did those with neocortical'...
    ' epilepsy (Supplemental Fig. 1B, Supplemental Table 3). There was no clear preictal change in spikes (Supplemental Results, Fig. 3A).</p>']);

%% Pre-ictal
all_pts_spikes_bins = sz_out.all_pts_spikes_bins;

nbins = size(all_pts_spikes_bins,2);
early_pre = 1:nbins/4; % first quarter
late_pre = nbins/4+1:nbins/2; % second quarter
early_late = [nanmean(all_pts_spikes_bins(:,early_pre),2),nanmean(all_pts_spikes_bins(:,late_pre),2)];


% Early vs late preictal spike rates
stats = plot_paired_data(early_late',{'early preictal','late preictal','late'},'Spikes/elec/min','paired',plot_type,0,0);

% Results text
fprintf(sfid,'<p><b>Preictal spike changes</b><br>');
fprintf(sfid,['<p>Visually, there was no obvious preictal increase in spikes (Fig 3A).'...
    ' Across all patients, the median spike rate in the late preictal period (median %1.2f spikes/elecs/min)'...
    ' was similar to that of the early preictal period (median %1.2f spikes/elecs/min) when examining 6-hour preictal periods '...
    '(Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Supplemental Fig 3C).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

% Secondary analyses
fprintf(sfid,[' However, spike rates were significantly higher in the late preictal period than'...
    ' in the early preictal period when performing our secondary analysis of 3-hour preictal windows (<i>p</i> = 0.018) '...
    'and 12-hour preictal windows (<i>p</i> = 0.017) (Supplemental Table 2).']); 

%% Rate diff in preictal period by epilepsy localization
rate_diff = early_late(:,2)-early_late(:,1);
[p,~,stats] = ranksum(rate_diff(temporal),rate_diff(extra));
% text
W = stats.ranksum;
nt = sum(~(isnan(rate_diff(temporal))));
ne = sum(~(isnan(rate_diff(extra))));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);

fprintf(sfid,[' For the primary analysis of a 6-hour preictal period, the early-to-late'...
    ' preictal spike rate change was similar between patients with temporal lobe epilepsy '...
    '(median = %1.2f spikes/elecs/min) and those with extra-temporal lobe epilepsy'...
    ' (median = %1.2f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>temporal</sub></i> = %d, <i>N<sub>extra-temporal</sub></i> = %d) ='...
    ' %1.1f, %s) (Supplemental Fig 3D). There was also no significant difference in the preictal spike '...
    'rate change between epilepsy localizations when studying the 3-hour '...
    '(<i>p</i> = 0.64) or 12-hour preictal windows (<i>p</i> = 0.41) (Supplemental Table 2).'],nanmedian(rate_diff(temporal)),nanmedian(rate_diff(extra)),...
    nt,ne,U,get_p_html(p));



%% Percent detected asleep
nexttile([1 1])
plot(times,mean_sleep,'color',[0.9290, 0.6940, 0.1250],'linewidth',2)
hold on

set(gca,'fontsize',15)
ylabel('% Classified asleep')
%title('Peri-ictal sleep classification')
xlabel('Hours surrounding seizure')
plot([0 0],ylim,'k--','linewidth',3)

%% Add annotations
annotation('textbox',[0 0.915 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.34 0.915 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0.685 0.915 0.1 0.1],'String','C','fontsize',25,'linestyle','none')


%% Add box around D
%annotation('rectangle',[0.001 0.001 0.32 0.65],'color','k','linewidth',2,'linestyle','--')

print([out_folder,'Fig3'],'-depsc')
close(gcf)
fclose(fid)
fclose(sfid)

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