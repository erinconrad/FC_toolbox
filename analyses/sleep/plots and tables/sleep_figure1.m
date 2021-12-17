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
set(gcf,'position',[100 100 1300 1100])
tiledlayout(3,3,'tilespacing','tight','padding','tight')

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

%% 2A overall spike rates
all_rate = bin_out.all_rates;
nexttile
plot_paired_data(all_rate(:,1:2)',{'wake','sleep'},'Spikes/elec/min','paired',plot_type)
title('Overall spike rate')

%% 2B independent spikes
seq_sw = bin_out.seq_sw;
nexttile
plot_paired_data(seq_sw(:,1:2)',{'wake','sleep'},'Spikes/elec/min','paired',plot_type)
title('Independent spikes')

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
title('Sleep-wake spike spread pattern correlation')
set(gca,'fontsize',15)
xl = xlim;
yl = ylim;
text(xl(2),yl(1),sprintf('Median r = %1.2f',nanmedian(rl_sw_corr)),...
    'fontsize',15,'verticalalignment','bottom','horizontalalignment','right')

%% 2E NS
ns_sw = bin_out.ns_sw;
nexttile
plot_paired_data(ns_sw',{'wake','sleep'},'Average node strength','paired',plot_type)
title('Functional connectivity')

print([out_folder,'Fig1'],'-dpng')

end