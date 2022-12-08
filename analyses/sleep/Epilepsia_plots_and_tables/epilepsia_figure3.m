function epilepsia_figure3(null_flag)

%% Parameters
min_rate = 0.1;
plot_type = 'scatter';
fsize = 20;

%{
myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];
%}
myColours = [0 33 87;...
122 80 113;...    
227 124 29;...
    86 152 163]/255;


locations = fc_toolbox_locs;
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Make fake data
if null_flag
    out = generate_fake_null_data(out);
end

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/epilepsia/'];

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');

figure
set(gcf,'position',[100 100 1100 650])
tiledlayout(2,3,'tilespacing','tight','padding','tight')

%% 2B independent spikes
seq_sw = bin_out.seq_sw; %1,2,3,4 are #seq wake, #seq sleep, seq length wake, seq length sleep, respectively
nexttile
stats = plot_paired_data(seq_sw(:,1:2)',{'wake','sleep'},'# Spike sequences','paired',plot_type);
title('Spike sequences')

% Results text
fprintf(fid,['<p>We next asked whether the higher spike rate in sleep was '...
    'driven by an increase in the number of spike sequences or a'...
    ' greater spread of each spike sequence (each spike sequence involving more electrodes).']);
fprintf(fid,[' The number of spike sequences was higher in sleep (median %1.1f spike sequences/min)'...
    ' than wake (median %1.1f spike sequences/min) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s, effect size r = %1.2f) (Fig. 3A).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval),stats.r);

%% 2C spike spread
nexttile
stats = plot_paired_data(seq_sw(:,3:4)',{'wake','sleep'},'# Spikes/sequence','paired',plot_type);
title('Spike spread')

% Results text
fprintf(fid,[' The spike spread was also higher in sleep (median %1.1f spikes/sequence)'...
    ' than wake (median %1.1f spikes/sequence) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s, effect size r = %1.2f) (Fig. 3B).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval),stats.r);

%{
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

nexttile
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
    'in spike spread patterns to be the Spearman correlation between the spike spread vectors '...
    'in sleep and that in wake. The spike rates and spread patterns '...
    'were highly consistent between'...
    ' wake and sleep (median &#961 across patients: %1.2f for rate, %1.2f for spread) (Fig. 3C).'],...
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
    %}

%% 2E NS
ns_sw = bin_out.ns_sw;
nexttile
stats = plot_paired_data(ns_sw',{'wake','sleep'},'Average node strength','paired',plot_type);
title('Functional connectivity')
fprintf(fid,[' To test a potential mechanism for the increased spike rates in sleep,'...
    ' we compared functional connectivity between wake and sleep. The functional connectivity'...
    ' as measured by the average node strength was higher in sleep (median %1.1f)'...
    ' than wake (median %1.1f) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s, effect size r = %1.2f) (Fig. 3C). '],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval),stats.r);

%% Prep sleep state stuff
all_ss = out.seeg_out.all_ss;
fc_ss = out.seeg_out.fc_ss;
spread_ss = out.seeg_out.spread_ss;
nseq_ss = out.seeg_out.nseq_ss;

% rm n1
n1 = strcmp(all_ss,'N1');
fc_ss(:,n1) = [];
spread_ss(:,n1) = [];
nseq_ss(:,n1) = [];
all_ss(n1) = [];

%% Spikes sequences by sleep state
nexttile
rows_with_any_nans_seq = any(isnan(nseq_ss),2);
[p,tbl_nseq] = friedman(nseq_ss(~rows_with_any_nans_seq,:),1,'off');
b=boxplot(nseq_ss,'labels',all_ss,'colors',myColours);

hold on
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
new_y = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 4],[ybar ybar],'k-','linewidth',2)
text(2.5,ytext,get_asterisks(p,1),'horizontalalignment','center','fontsize',20)
ylim(new_y)

set(b,'linewidth',2)
set(gca,'fontsize',fsize)
ylabel('# Spike sequences')
title('Spike sequences across stages')

h = findobj(gcf,'tag','Outliers');
for i = 1:numel(h)
    xpos = h(i).XData(1);
    h(i).MarkerEdgeColor = myColours(xpos,:);
end

%% Spikes spread by sleep state
nexttile
rows_with_any_nans_spread = any(isnan(spread_ss),2);
[p,tbl_spread] = friedman(spread_ss(~rows_with_any_nans_seq,:),1,'off');
b=boxplot(spread_ss,'labels',all_ss,'colors',myColours);

hold on
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
new_y = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 4],[ybar ybar],'k-','linewidth',2)
text(2.5,ytext,get_asterisks(p,1),'horizontalalignment','center','fontsize',20)
ylim(new_y)

set(b,'linewidth',2)
set(gca,'fontsize',fsize)
ylabel('# Spikes/sequence')
title('Spike spread across stages')

h = findobj(gcf,'tag','Outliers');
for i = 1:numel(h)
    xpos = h(i).XData(1);
    h(i).MarkerEdgeColor = myColours(xpos,:);
end

%% FC by sleep state
nexttile
rows_with_any_nans_fc = any(isnan(fc_ss),2);
[p,tbl_fc] = friedman(fc_ss(~rows_with_any_nans_seq,:),1,'off');
b=boxplot(fc_ss,'labels',all_ss,'colors',myColours);
set(b,'linewidth',2)
set(gca,'fontsize',fsize)
ylabel('Average node strength')
title('Connectivity across stages')


hold on
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
new_y = [yl(1) yl(1) + 1.2*(yl(2)-yl(1))];
plot([1 4],[ybar ybar],'k-','linewidth',2)
text(2.5,ytext,get_asterisks(p,1),'horizontalalignment','center','fontsize',20)
ylim(new_y)

h = findobj(gcf,'tag','Outliers');
for i = 1:numel(h)
    xpos = h(i).XData(1);
    if isnan(xpos), continue; end
    h(i).MarkerEdgeColor = myColours(xpos,:);
end



fprintf(fid,['The number of spike sequences, spike spread, and functional connectivity '...
    'all varied across specific sleep stages, again demonstrating generally higher '...
    'values in NREM sleep (Friedman test: p < 0.001 for each analysis; '...
    'Figure 3D-F; Tables S3-5 show descriptive statistics and pairwise comparisons). <p>'])

sleep_comparison_table(nseq_ss,all_ss,'TableS3.html')
sleep_comparison_table(spread_ss,all_ss,'TableS4.html')
sleep_comparison_table(fc_ss,all_ss,'TableS5.html')


%% Add annotations
annotation('textbox',[0 0.905 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.33 0.905 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0.655 0.905 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.405 0.1 0.1],'String','D','fontsize',25,'linestyle','none')
annotation('textbox',[0.33 0.405 0.1 0.1],'String','E','fontsize',25,'linestyle','none')
annotation('textbox',[0.655 0.405 0.1 0.1],'String','F','fontsize',25,'linestyle','none')
fontname(gcf,"calibri");
fclose(fid);
print([out_folder,'Fig3'],'-depsc')
print([out_folder,'Fig3'],'-dpng')
close(gcf)


end