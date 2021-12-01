function assess_soz_spike_ranking

%% Parameters
nb = 1e4;

locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
out = load([out_folder,'out.mat']);
out = out.out;

%% Get info
rates = out.bin_out.all_elec_rates;
soz = out.bin_out.all_is_soz;
locs = out.circ_out.all_locs;
is_temporal = cellfun(@(x) strcmp(x,'temporal'),locs);

% convert soz to elec nums rather than bin
soz = cellfun(@find,soz,'uniformoutput',false);

%% Do MC test
mcout = test_ranking(rates,soz,nb);
median_ranking_mc = mcout.median_ranking_mc;
median_ranking_true = mcout.median_ranking_true;
all_rankings = mcout.all_rankings;
pval = mcout.pval;


%% Plots
figure
set(gcf,'position',[10 10 1400 400])
tiledlayout(1,3,'padding','compact','tilespacing','compact')

nexttile
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
ylabel('Median spike rate ranking')
title('Spike rate ranking of SOZ compared to chance')

nexttile
plot(1+randn(sum(is_temporal),1)*0.05,all_rankings(is_temporal),'o','linewidth',2)
hold on
plot(2+randn(sum(~is_temporal),1)*0.05,all_rankings(~is_temporal),'o','linewidth',2)
pval = ranksum(all_rankings(is_temporal),all_rankings(~is_temporal));
yl = ylim;
ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2)-yl(1));
newyl = [yl(1) yl(1) + 1.17*(yl(2)-yl(1))];
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,ytext,get_p_text(pval),'horizontalalignment','center','fontsize',15)
xticks([1 2])
xticklabels({'Temporal','Extra-temporal'})
set(gca,'fontsize',15)
ylabel('Median SOZ spike rate ranking')
title('Spike rate ranking of SOZ by seizure localization')

nexttile
soz_rank_sw = out.bin_out.soz_rank_sw;
plot_paired_data(soz_rank_sw',{'wake','sleep'},'Rank in spike rate','paired','scatter','ranking')
title({'SOZ spike rate ranking','by wake vs sleep'})

print([out_folder,'Fig4'],'-dpng')

end