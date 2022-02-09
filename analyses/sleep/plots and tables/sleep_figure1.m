function sleep_figure1

%% Parameters
min_rate = 0.1;
plot_type = 'scatter';
nblocks = 6;
%{
myColours = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250];
%}

myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];


locations = fc_toolbox_locs;
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');

fprintf(fid,'<p><b>Changes in spikes with sleep</b><br>');

figure
set(gcf,'position',[100 100 1200 1100])
tiledlayout(3,6,'tilespacing','compact','padding','compact')

%% A - PSD
median_psd = circ_out.median_psd;
iqr_psd = circ_out.iqr_psd;
periods = circ_out.periods;

nexttile([1 2])
shaded_error_bars_fc(periods,median_psd,iqr_psd,myColours(1,:));
xlim([0 100])
xlabel('Period (hours)')
ylabel({'Spike rate power index'});
set(gca,'fontsize',15)
title('Spike rate periodogram')

fprintf(fid,['We asked how spike rates changed with sleep.'...
    ' Our analysis of the circadian power of spike rates demonstrated a '...
    'clear visual peak at 24 hours, consistent with the hypothesis that '...
    'spike rates demonstrate a circadian rhythm (Fig. 2A).']); 

%% TOD
all_tod_rate = circ_out.all_tod_rate;
tod_edges = bin_out.tod_edges;
counts = nansum(all_tod_rate,1);
hours_mins = convert_edges_to_hours_mins(tod_edges);
nbins = length(counts);
skip = 18;
norm_rate = (all_tod_rate - nanmedian(all_tod_rate,2))./iqr(all_tod_rate,2);
nexttile([1 2])
median_tod_rate = (nanmedian(norm_rate,1));
polar = convert_times_to_polar(tod_edges,'radians');
%
all_tod_sw = out.bin_out.all_tod_sw;
ind_pt_prop = all_tod_sw(:,:,2)./(all_tod_sw(:,:,2)+all_tod_sw(:,:,1));
prop_asleep = squeeze(nanmean(ind_pt_prop,1))*3.5;
psleep = polarhistogram('BinEdges',polar,'BinCounts',prop_asleep,...
    'displayStyle','stairs','linewidth',1,'edgecolor',[0.9290, 0.6940, 0.1250],'linestyle','-');
hold on
% Plot median of spike rate over time
mspikes = repmat(nanmedian(median_tod_rate+min(median_tod_rate)+1),1,length(median_tod_rate));
%
pm = polarhistogram('BinEdges',polar,'BinCounts',mspikes,...
    'displayStyle','stairs','linewidth',1,'edgecolor',myColours(1,:),'linestyle','--');
pspikes = polarhistogram('BinEdges',polar,'BinCounts',median_tod_rate+min(median_tod_rate)+1,...
    'displayStyle','stairs','linewidth',2,'edgecolor',myColours(1,:));


%}
%}
%{
polarhistogram('BinEdges',polar,'BinCounts',counts,...
    'edgecolor','none')
%}
% Add in sw



set(gca,'ThetaDir','clockwise');
set(gca,'ThetaZeroLocation','top');
set(gca,'rticklabels',[])
thetaticks(polar(1:skip:nbins)*360/(2*pi))
thetaticklabels(hours_mins(1:skip:nbins+1))
set(gca,'fontsize',15)
title('Normalized spike rate and % asleep')
lp = legend([pspikes,pm,psleep],{'Spike rates','Median spike rate','% asleep'},'fontsize',15,...
    'Position',[0.3360 0.8300 0.1333 0.0550],'box','off');


observations = convert_counts_to_observations(counts,tod_edges);
polar2 = convert_times_to_polar(observations,'radians');
%circ_plot(polar,'hist',[],length(tod_edges),true,true,'linewidth',2,'color','r')
%[pval z] = circ_rtest(polar2);


[pval z all_mu] = test_pt_circular_means(all_tod_rate,polar,hours_mins);
%circ_plot(all_mu,'hist',[],length(polar),true,true,'linewidth',2,'color','r')


tl = thetalim;
rl = rlim;
text(pi,1,get_p_text(pval),'fontsize',15,'horizontalalignment','center');

fprintf(fid,[' Examination of spike rates by time of day revealed a non-uniform distribution'...
    ' (Rayleigh test of patients'' circular means: z = %1.1f, %s). Visual analysis demonstrated '...
    'that spike rates tend to peak in the early morning and have a nadir in the late morning,'...
    ' afternoon, and evening (Fig. 2B).'],z,get_p_html(pval)); 

%% 2A overall spike rates
all_rate = bin_out.all_rates;
nexttile([1 2])
stats = plot_paired_data(all_rate(:,1:2)',{'wake','sleep'},'Spikes/elec/min','paired',plot_type);
title('Spike rate in wake and sleep')

% Results text
fprintf(fid,[' Across patients, spike rates were higher in sleep (median %1.1f spikes/elecs/min)'...
    ' than wake (median %1.1f spikes/elecs/min) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s).'...
    ' Spike rates were higher in sleep for %d out of %d patients (Fig. 2C).</p>'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval),...
    stats.nhigher_n(1),stats.nhigher_n(2));

%% 2B independent spikes
seq_sw = bin_out.seq_sw; %1,2,3,4 are #seq wake, #seq sleep, seq length wake, seq length sleep, respectively
nexttile([1 2])
stats = plot_paired_data(seq_sw(:,1:2)',{'wake','sleep'},'# Spike sequences','paired',plot_type);
title('Independent spikes')

% Results text
fprintf(fid,['<p>We next asked whether the higher spike rate in sleep was '...
    'driven by an increase in the number of independent spike sequences or a'...
    ' greater spread of individual spike sequences (each spike sequence involving more electrodes).']);
fprintf(fid,[' The number of independent spike sequences was higher in sleep (median %1.1f spike sequences/min)'...
    ' than wake (median %1.1f spike sequences/min) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 2D).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

%% 2C spike spread
nexttile([1 2])
stats = plot_paired_data(seq_sw(:,3:4)',{'wake','sleep'},'# Spikes/sequence','paired',plot_type);
title('Spike spread')

% Results text
fprintf(fid,[' The spike spread was also higher in sleep (median %1.1f spikes/sequence)'...
    ' than wake (median %1.1f spikes/sequence) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 2E).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

%% 2D RL sleep-wake correlation across patients
%rl_sw_corr = bin_out.rl_sw_corr;

% Get electrodes above minimum rate (for purpose of RL analysis)
all_elecs_rates = bin_out.all_elecs_rates;
npts = length(all_elecs_rates);
min_rate_elecs = cell(npts,1);
for i = 1:npts
    curr = all_elecs_rates{i};
    min_rate_elecs{i} = curr >= min_rate;
end

% SW consistency in spike rates (spearman correlation between rates in wake
% and rates in sleep)
all_elecs_rates_sw = bin_out.all_elecs_rates_sw;
rates_sw_corr = cellfun(@(x) corr(x(:,1),x(:,2),'type','spearman','rows','pairwise'),...
    all_elecs_rates_sw);

% Sw consistency in spike spread (remove those with few spikes)
all_elecs_rl_sw = bin_out.all_elecs_rl_sw;
rl_sw_corr = nan(npts,1);
for i = 1:npts
    curr = all_elecs_rl_sw{i};
    curr_min = min_rate_elecs{i};
    curr = curr(curr_min,:);
    rl_sw_corr(i) = corr(curr(:,1),curr(:,2),'type','spearman','rows','pairwise');
end

nexttile([1 2])
plot(1+randn(npts,1)*0.05,rates_sw_corr,'o','linewidth',2,'color',myColours(1,:))
hold on
plot([0.75 1.25],[nanmedian(rates_sw_corr) nanmedian(rates_sw_corr)],'linewidth',2,'color',myColours(1,:))
plot(2+randn(npts,1)*0.05,rl_sw_corr,'o','linewidth',2,'color',[0.6350, 0.0780, 0.1840])
plot([1.75 2.25],[nanmedian(rl_sw_corr) nanmedian(rl_sw_corr)],'linewidth',2,'color',[0.6350, 0.0780, 0.1840])
ylabel('Correlation coefficient')
xlim([0.5 2.5])
plot(xlim,[0 0],'k--','linewidth',2)
xticks([1 2])
xticklabels({'Spike rate','Spike spread'})
ylim([-1 1])
title('Consistency in spikes from wake-to-sleep')
set(gca,'fontsize',15)

fprintf(fid,[' We measured the consistency in the patterns of spike rates and spike spread '...
    'from sleep to wake. We defined the sleep-wake consistency in spike rate '...
    'patterns to be the Spearman correlation between the vectors of spike rates '...
    'across electrodes for wake and sleep. We defined the sleep-wake consistency '...
    'in spike spread patterns to be the Spearman correlation between the spike timing vectors '...
    'in sleep and that in wake. In the timing analysis, we excluded electrodes with average spike'...
    ' rates less than 0.1 spikes/minute in order to prevent artifactual detections from'...
    ' contributing to the timing analysis. The spike rates and spread patterns '...
    'were highly consistent between'...
    ' wake and sleep (median &#961 across patients: %1.2f for rate, %1.2f for spread) (Fig. 2F).'],...
    nanmedian(rates_sw_corr),nanmedian(rl_sw_corr));
%{
nexttile([1 2])
plot(rl_sw_corr,'o','linewidth',2)
hold on
plot(xlim,[0 0],'k--','linewidth',2)
ylabel('Correlation coefficient')
xlabel('Patient')
xticklabels([])
ylim([-1 1])
title('Sleep-wake spike spread consistency')
set(gca,'fontsize',15)
xl = xlim;
yl = ylim'
tV"
t;
text(xl(2),yl(1),sprintf('Median r = %1.2f',nanmedian(rl_sw_corr)),...
    'fontsize',15,'verticalalignment','bottom','horizontalalignment','right')
%}

%% 2E NS
ns_sw = bin_out.ns_sw;
nexttile([1 2])
stats = plot_paired_data(ns_sw',{'wake','sleep'},'Average node strength','paired',plot_type);
title('Functional connectivity')
fprintf(fid,[' To test a potential mechanism for the increased spike rates in sleep,'...
    ' we compared functional connectivity between wake and sleep. The functional connectivity'...
    ' as measured by the average node strength was higher in sleep (median %1.1f)'...
    ' than wake (median %1.1f) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 2G).</p>'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

%% Overall rate localization
rate = bin_out.overall_rates;
loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');

if 0
   table(loc,temporal,extra) 
end

[p,~,stats] = ranksum(rate(temporal),rate(extra));

nexttile([1 2])
plot(1+randn(sum(temporal),1)*0.05,rate(temporal),'o','linewidth',2,'color',myColours(2,:))
hold on
plot([0.7 1.3],[nanmedian(rate(temporal)) nanmedian(rate(temporal))],...
    'linewidth',2,'color',myColours(2,:))
plot(2+randn(sum(extra),1)*0.05,rate(extra),'o','linewidth',2,'color',myColours(3,:))
plot([1.7 2.3],[nanmedian(rate(extra)) nanmedian(rate(extra))],...
    'linewidth',2,'color',myColours(3,:))
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})
ylabel('Spikes/elec/min')
title('Average spike rate by localization')
set(gca,'fontsize',15);
xlim([0 3])
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.13*(yl(2)-yl(1));
ylnew = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
ylim(ylnew)

W = stats.ranksum;
nt = sum(~isnan(rate(temporal)));
ne = sum(~isnan(rate(extra)));
ns = min([sum(temporal),sum(extra)]);
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,['<p>We next tested how overall spike rates varied by epilepsy localization. There was no difference in overall spike rates between patients'...
    ' with temporal lobe epilepsy (median = %1.1f spikes/elecs/min) and patients with extra-temporal'...
    ' lobe epilepsy (median = %1.1f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>temporal</sub></i> = %d, <i>N<sub>extra-temporal</sub></i> = %d) ='...
    ' %1.1f, %s) (Fig. 2H).'],nanmedian(rate(temporal)),nanmedian(rate(extra)),...
    nt,ne,U,get_p_html(p));

%% Add Supplemental Figure 1 result
mt = contains(loc,'mesial temporal');
nc = contains(loc,'cortical') | contains(loc,'other cortex'); % temporal neocortical and other cortex
[p,~,stats] = ranksum(rate(mt),rate(nc));
W = stats.ranksum;
nt = sum(~isnan(rate(mt)));
ne = sum(~isnan(rate(nc)));
ns = min([sum(mt),sum(nc)]);
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,[' There was also no difference in overall spike rates between patients'...
    ' with mesial temporal lobe epilepsy (median = %1.1f spikes/elecs/min) and patients with neocortical'...
    ' epilepsy (median = %1.1f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>mesial temporal</sub></i> = %d, <i>N<sub>neocortical</sub></i> = %d) ='...
    ' %1.1f, %s) (Supplementary Fig. 1A).</p>'],nanmedian(rate(mt)),nanmedian(rate(nc)),...
    nt,ne,U,get_p_html(p));


%% 2F Localization
rate_sw = bin_out.all_rates;
loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
rate_diff = (rate_sw(:,2)-rate_sw(:,1));
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
[p,~,stats] = ranksum(rate_diff(temporal),rate_diff(extra));

W = stats.ranksum;
nt = sum(~isnan(rate_diff(temporal)));
ne = sum(~isnan(rate_diff(extra)));
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);
fprintf(fid,['<p>Finally, we tested whether the sleep-wake spike rate difference varied by epilepsy localization. '...
    'There was a non-significant trend toward higher sleep-wake spike rate increase in patients'...
    ' with temporal lobe epilepsy (median = %1.1f spikes/elecs/min) compared to patients with extra-temporal'...
    ' lobe epilepsy (median = %1.1f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>temporal</sub></i> = %d, <i>N<sub>extra-temporal</sub></i> = %d) ='...
    ' %1.1f, %s) (Fig. 2I).'],nanmedian(rate_diff(temporal)),nanmedian(rate_diff(extra)),...
    nt,ne,U,get_p_html(p));
%{
mesial_temporal = strcmp(loc,'mesial temporal');
neocortical = strcmp(loc,'temporal neocortical') | strcmp(loc,'other cortex');
diffuse = strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
nmt = sum(mesial_temporal);
nn = sum(neocortical);
nd = sum(diffuse);

plot(1+randn(nmt,1)*0.05,rate_diff(mesial_temporal),'o')
hold on
plot(2+randn(nn,1)*0.05,rate_diff(neocortical),'o')
plot(3+randn(nd,1)*0.05,rate_diff(diffuse),'o')
%}

nexttile([1 2])
plot(1+randn(sum(temporal),1)*0.05,rate_diff(temporal),'o','linewidth',2,'color',myColours(2,:))
hold on
plot([0.7 1.3],[nanmedian(rate_diff(temporal)) nanmedian(rate_diff(temporal))],...
    'linewidth',2,'color',myColours(2,:))
plot(2+randn(sum(extra),1)*0.05,rate_diff(extra),'o','linewidth',2,'color',myColours(3,:))
plot([1.7 2.3],[nanmedian(rate_diff(extra)) nanmedian(rate_diff(extra))],...
    'linewidth',2,'color',myColours(3,:))
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

%% Add Supplemental result
[p,~,stats] = ranksum(rate_diff(mt),rate_diff(nc));
W = stats.ranksum;
nt = sum(~isnan(rate_diff(mt)));
ne = sum(~isnan(rate_diff(nc)));
ns = min([sum(mt),sum(nc)]);
U1 = W - nt*(nt+1)/2;
U2 = nt*ne-U1;
U = min([U1,U2]);

fprintf(fid,[' There was no significant difference in sleep-wake spike rate difference between patients'...
    ' with mesial temporal lobe epilepsy (median = %1.1f spikes/elecs/min) and patients with neocortical'...
    ' epilepsy (median = %1.1f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>mesial temporal</sub></i> = %d, <i>N<sub>neocortical</sub></i> = %d) ='...
    ' %1.1f, %s) (Supplementary Fig. 1B).</p>'],nanmedian(rate_diff(mt)),nanmedian(rate_diff(nc)),...
    nt,ne,U,get_p_html(p));

%% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.33 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0.65 0.91 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.58 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0.33 0.58 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0.65 0.58 0.1 0.1],'String','F','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.24 0.1 0.1],'String','G','fontsize',25,'linestyle','none')
annotation('textbox',[0.33 0.24 0.1 0.1],'String','H','fontsize',25,'linestyle','none')
annotation('textbox',[0.65 0.24 0.1 0.1],'String','I','fontsize',25,'linestyle','none')

fclose(fid);
print([out_folder,'Fig2'],'-depsc')

end