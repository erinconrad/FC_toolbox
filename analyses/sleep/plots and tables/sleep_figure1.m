function sleep_figure1

%% Parameters
plot_type = 'scatter';
nblocks = 6;
myColours = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250];



locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

%% Load out file and get roc stuff
out = load([out_folder,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];


figure
set(gcf,'position',[100 100 1100 1100])
tiledlayout(3,6,'tilespacing','compact','padding','compact')

%% A - PSD
median_psd = circ_out.median_psd;
iqr_psd = circ_out.iqr_psd;
periods = circ_out.periods;

nexttile([1 3])
shaded_error_bars(periods,median_psd,iqr_psd,[]);
xlim([0 100])
xlabel('Period (hours)')
ylabel({'Spike rate power index'});
set(gca,'fontsize',15)
title('Spike rate periodogram')

%% 2A overall spike rates
all_rate = bin_out.all_rates;
nexttile([1 3])
plot_paired_data(all_rate(:,1:2)',{'wake','sleep'},'Spikes/elec/min','paired',plot_type)
title('Spike rate in wake and sleep')

%% 2B independent spikes
seq_sw = bin_out.seq_sw;
nexttile([1 2])
plot_paired_data(seq_sw(:,1:2)',{'wake','sleep'},'# Spike sequences','paired',plot_type)
title('Independent spikes')

%% 2C spike spread
nexttile([1 2])
plot_paired_data(seq_sw(:,3:4)',{'wake','sleep'},'# Spikes/sequence','paired',plot_type)
title('Spike spread')

%% 2D RL sleep-wake correlation across patients
rl_sw_corr = bin_out.rl_sw_corr;
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
yl = ylim;
text(xl(2),yl(1),sprintf('Median r = %1.2f',nanmedian(rl_sw_corr)),...
    'fontsize',15,'verticalalignment','bottom','horizontalalignment','right')

%% 2E NS
ns_sw = bin_out.ns_sw;
nexttile([1 2])
plot_paired_data(ns_sw',{'wake','sleep'},'Average node strength','paired',plot_type)
title('Functional connectivity')

%% Overall rate localization
rate = bin_out.overall_rates;
loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
p = ranksum(rate(temporal),rate(extra));

nexttile([1 2])
plot(1+randn(sum(temporal),1)*0.05,rate(temporal),'o','linewidth',2,'color',myColours(1,:))
hold on
plot([0.7 1.3],[nanmedian(rate(temporal)) nanmedian(rate(temporal))],...
    'linewidth',2,'color',myColours(1,:))
plot(2+randn(sum(extra),1)*0.05,rate(extra),'o','linewidth',2,'color',myColours(2,:))
plot([1.7 2.3],[nanmedian(rate(extra)) nanmedian(rate(extra))],...
    'linewidth',2,'color',myColours(2,:))
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


%{
% Correlate increase in NS with increase in spike rate
ns_inc = (ns_sw(:,2)-ns_sw(:,1))./ns_sw(:,1);
rate_inc = (all_rate(:,2)-all_rate(:,1))./all_rate(:,1);
[r,p] = corr(ns_inc,rate_inc,'type','spearman');
plot(ns_inc,rate_inc,'o','linewidth',2)
xlabel('Relative node strength increase in sleep')
ylabel('Relative spike rate increase in sleep')
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('\\rho = %1.2f, %s',r,get_p_text(p)),...
'verticalalignment','top','fontsize',15)
set(gca,'fontsize',15)

% Correlate electrode-specific increase in NS with increase in spike rate
elecs_ns = bin_out.all_elecs_ns_sw;
elecs_rate = bin_out.all_elecs_rates_sw;
npts = length(elecs_ns);
all_corr = nan(npts,1);
for i = 1:npts
    ns_change = elecs_ns{i}(:,2) - elecs_ns{i}(:,1);
    rate_change = elecs_rate{i}(:,2) - elecs_rate{i}(:,1);
    all_corr(i) = corr(ns_change,rate_change,'type','spearman','rows','pairwise');
end
plot(all_corr,'o')

%}

%% 2F Localization
rate_sw = bin_out.all_rates;
loc = circ_out.all_locs;
temporal = contains(loc,'temporal');
rate_diff = (rate_sw(:,2)-rate_sw(:,1));
extra = strcmp(loc,'other cortex') | strcmp(loc,'diffuse') | strcmp(loc,'multifocal');
p = ranksum(rate_diff(temporal),rate_diff(extra));

nexttile([1 2])
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

%% Add annotations
annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.48 0.91 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.58 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.32 0.58 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0.64 0.58 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.24 0.1 0.1],'String','F','fontsize',25,'linestyle','none')
annotation('textbox',[0.32 0.24 0.1 0.1],'String','G','fontsize',25,'linestyle','none')
annotation('textbox',[0.64 0.24 0.1 0.1],'String','H','fontsize',25,'linestyle','none')


print([out_folder,'Fig1'],'-dpng')

end