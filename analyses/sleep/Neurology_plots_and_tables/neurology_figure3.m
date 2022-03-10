function neurology_figure3

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
out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/neurology/'];

%% Prep output text file
fid = fopen([out_folder,'results.html'],'a');

figure
set(gcf,'position',[100 100 900 600])
tiledlayout(2,2,'tilespacing','compact','padding','compact')

%% 2B independent spikes
seq_sw = bin_out.seq_sw; %1,2,3,4 are #seq wake, #seq sleep, seq length wake, seq length sleep, respectively
nexttile
stats = plot_paired_data(seq_sw(:,1:2)',{'wake','sleep'},'# Spike sequences','paired',plot_type);
title('Independent spikes')

% Results text
fprintf(fid,['<p>We next asked whether the higher spike rate in sleep was '...
    'driven by an increase in the number of independent spike sequences or a'...
    ' greater spread of individual spike sequences (each spike sequence involving more electrodes).']);
fprintf(fid,[' The number of independent spike sequences was higher in sleep (median %1.1f spike sequences/min)'...
    ' than wake (median %1.1f spike sequences/min) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 3A).'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));

%% 2C spike spread
nexttile
stats = plot_paired_data(seq_sw(:,3:4)',{'wake','sleep'},'# Spikes/sequence','paired',plot_type);
title('Spike spread')

% Results text
fprintf(fid,[' The spike spread was also higher in sleep (median %1.1f spikes/sequence)'...
    ' than wake (median %1.1f spikes/sequence) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 3B).'],...
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

%% 2E NS
ns_sw = bin_out.ns_sw;
nexttile
stats = plot_paired_data(ns_sw',{'wake','sleep'},'Average node strength','paired',plot_type);
title('Functional connectivity')
fprintf(fid,[' To test a potential mechanism for the increased spike rates in sleep,'...
    ' we compared functional connectivity between wake and sleep. The functional connectivity'...
    ' as measured by the average node strength was higher in sleep (median %1.1f)'...
    ' than wake (median %1.1f) (Wilcoxon signed-rank test: <i>T<sup>+</sup></i> = %1.1f, %s) (Fig. 3D).</p>'],...
    stats.medians(2),stats.medians(1),stats.Tpos,get_p_html(stats.pval));


%% Add annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.9 0.1 0.1],'String','B','fontsize',25,'linestyle','none')
annotation('textbox',[0 0.41 0.1 0.1],'String','C','fontsize',25,'linestyle','none')
annotation('textbox',[0.5 0.41 0.1 0.1],'String','D','fontsize',25,'linestyle','none')

fclose(fid);
print([out_folder,'Fig3'],'-depsc')
close(gcf)


end