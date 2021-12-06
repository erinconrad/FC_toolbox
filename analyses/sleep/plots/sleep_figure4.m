function sleep_figure4

%% Seed rng (because it's a MC test)
rng(0)
nb = 1e4;
min_rate = 0.1;


%% Parameters
plot_type = 'scatter';
nblocks = 6;
myColours = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];



locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

%% Load out file and get roc stuff
out = load([out_folder,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

figure
set(gcf,'position',[100 100 900 700])
tiledlayout(2,6,'tilespacing','tight','padding','tight')

%% Get stuff
rate = out.bin_out.all_elec_rates;
rl = out.bin_out.all_elecs_rl;
soz_rank_sw_rate = out.bin_out.soz_rank_sw;
soz_rank_sw_rl = out.bin_out.soz_rank_sw_rl;
soz = out.bin_out.all_is_soz;
locs = out.circ_out.all_locs;
is_temporal = cellfun(@(x) strcmp(x,'temporal'),locs);

%% Do MC tests
mcout_rate = test_ranking(rate,soz,nb,'rate',rate,min_rate);
median_ranking_mc_rate = mcout_rate.median_ranking_mc;
median_ranking_true_rate = mcout_rate.median_ranking_true;
all_rankings_rate = mcout_rate.all_rankings;
pval_rate = mcout_rate.pval;

mcout_rl = test_ranking(rl,soz,nb,'rl',rate,min_rate);
median_ranking_mc_rl = mcout_rl.median_ranking_mc;
median_ranking_true_rl = mcout_rl.median_ranking_true;
all_rankings_rl = mcout_rl.all_rankings;
pval_rl = mcout_rl.pval;


%% SOZ spike rate ranking
thing_text = 'rate';
median_ranking_mc = median_ranking_mc_rate;
median_ranking_true = median_ranking_true_rate;
pval = pval_rate;
nexttile([1 2])
plot(sort(median_ranking_mc),'o','linewidth',2)
hold on
plot(xlim,[median_ranking_true median_ranking_true],'linewidth',2)
xl = xlim;
text(xl(1),median_ranking_true+2,...
    sprintf('Median ranking: %1.1f, %s',median_ranking_true,get_p_text(pval)),...
    'verticalalignment','bottom','fontsize',15)
legend({'Random electrodes','True SOZ electrodes'},'Location','Northwest','fontsize',15)
set(gca,'fontsize',15)
xticklabels([])
xlabel('Monte Carlo iteration')
ylabel(sprintf('Median spike %s ranking',thing_text))
title(sprintf('Spike %s ranking of SOZ compared to chance',thing_text))


%% SOZ spike timing ranking
thing_text = 'timing';
median_ranking_mc = median_ranking_mc_rl;
median_ranking_true = median_ranking_true_rl;
pval = pval_rl;
nexttile([1 2])
plot(sort(median_ranking_mc),'o','linewidth',2)
hold on
plot(xlim,[median_ranking_true median_ranking_true],'linewidth',2)
xl = xlim;
text(xl(1)+20,median_ranking_true+2,...
    sprintf('Median ranking: %1.1f, %s',median_ranking_true,get_p_text(pval)),...
    'verticalalignment','bottom','fontsize',15)
legend({'Random electrodes','True SOZ electrodes'},'Location','Northwest','fontsize',15)
set(gca,'fontsize',15)
xticklabels([])
xlabel('Monte Carlo iteration')
ylabel(sprintf('Median spike %s ranking',thing_text))
title(sprintf('Spike %s ranking of SOZ compared to chance',thing_text))

%% Correlate rate and rl
rl_rate_corr = nan(length(rate),1);
for ip = 1:length(rate)
    curr_rate = rate{ip};
    curr_rl = rl{ip};
    % Only do it for those with enough spikes
    spikey = curr_rate > min_rate;
    rl_rate_corr(ip) = corr(curr_rate(spikey),curr_rl(spikey),'type','spearman','rows','pairwise');
end

% Plot it
nexttile([1 2])
plot(rl_rate_corr,'o','linewidth',2)
hold on
plot(xlim,[0 0],'k--','linewidth',2)
xticklabels([])
xlabel('Patient')
ylabel('Correlation coefficients');
set(gca,'fontsize',15)
title({'Correlation between spike rate and spike timing'})

%% Spike rate ranking sleep vs wake
nexttile([1 3])
plot_paired_data(soz_rank_sw_rate',{'wake','sleep'},sprintf('Rank in spike rate'),'paired','scatter','ranking')
title(sprintf('SOZ spike rate ranking\nby wake vs sleep'))

%% Spike timing ranking sleep vs wake
nexttile([1 3])
plot_paired_data(soz_rank_sw_rl',{'wake','sleep'},sprintf('Rank in spike timing'),'paired','scatter','ranking')
title(sprintf('SOZ spike timing ranking\nby wake vs sleep'))

end