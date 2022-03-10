function neurology_figure2

%% Parameters
min_rate = 0.1;
plot_type = 'scatter';


myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];


locations = fc_toolbox_locs;
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/neurology/'];

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');

fprintf(fid,'<p><b>Changes in spikes with sleep</b><br>');

figure
set(gcf,'position',[10 10 900 600])
tiledlayout(2,2,'tilespacing','tight','padding','tight')

%% A - PSD
median_psd = circ_out.median_psd;
iqr_psd = circ_out.iqr_psd;
periods = circ_out.periods;

nexttile
shaded_error_bars_fc(periods,median_psd,iqr_psd,myColours(1,:));
xlim([0 100])
xlabel('Period (hours)')
ylabel({'Spike rate power index'});
set(gca,'fontsize',15)
title('Spike rate periodogram')

fprintf(fid,['Our analysis of the circadian power of spike rates demonstrated a '...
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
nexttile
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
    'Position',[0.52 0.8100 0.1333 0.0550],'box','off');


observations = convert_counts_to_observations(counts,tod_edges);
polar2 = convert_times_to_polar(observations,'radians');
%circ_plot(polar,'hist',[],length(tod_edges),true,true,'linewidth',2,'color','r')
%[pval z] = circ_rtest(polar2);


[pval,z,all_mu] = test_pt_circular_means(all_tod_rate,polar,hours_mins);
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
nexttile
stats = plot_paired_data(all_rate(:,1:2)',{'wake','sleep'},'Spikes/elec/min','paired',plot_type);
title('Spike rate in wake and sleep')

% Results text
fprintf(fid,[' Across patients, spike rates were higher in sleep (median %1.1f spikes/elecs/min)'...
    ' than wake (median %1.1f spikes/elecs/min) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s).'...
    ' (Fig. 2C).</p>'],...
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
    ' %1.1f, %s).'],nanmedian(rate(temporal)),nanmedian(rate(extra)),...
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
fprintf(fid,['We also tested whether the sleep-wake spike rate difference varied by epilepsy localization. '...
    'Patients with temporal lobe epilepsy had a higher sleep-wake spike rate increase '...
    '(median = %1.1f spikes/elecs/min) compared to patients with extra-temporal'...
    ' lobe epilepsy (median = %1.1f spikes/elecs/min) (Mann-Whitney test: <i>U</i>'...
    '(<i>N<sub>temporal</sub></i> = %d, <i>N<sub>extra-temporal</sub></i> = %d) ='...
    ' %1.1f, %s) (Fig. 2D).'],nanmedian(rate_diff(temporal)),nanmedian(rate_diff(extra)),...
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

nexttile
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


fprintf(fid,[' There was no significant difference in sleep-wake spike rate difference between patients'...
    ' with mesial temporal lobe epilepsy and patients with neocortical epilepsy '...
    '(Supplementary Fig. 1A, Supplementary Table 2).</p>']);

%% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.41 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.41 0.1 0.1],'String','D','fontsize',25,'linestyle','none')

fclose(fid);
print([out_folder,'Fig2'],'-depsc')
close(gcf)


end