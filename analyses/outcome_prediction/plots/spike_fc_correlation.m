function spike_fc_correlation

%% Parameters
do_r2 = 0;

locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
atlas_folder = [results_folder,'analysis/atlas/'];


%% Load out file and get roc stuff
out = load([out_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
rl = out.all_rl;
ns = out.all_ns;
soz = out.all_soz_bin;
soz_loc = out.all_soz_locs;
npts = length(soz);
labels = out.all_labels;
locs = out.all_locs;
fc = out.all_fc;

if do_r2
    fc = cellfun(@(x) x.^2, fc,'uniformoutput',false);
end

%% Get spike-fc correlations
spike_fc_corr = nan(npts,1);
for ip = 1:npts
    curr_labels = labels{ip};
    assert(sum(find_non_intracranial(curr_labels)) == 0)
    
    curr_avg_fc = nanmean(fc{ip},2);
    curr_spikes = rate{ip};
    spike_fc_corr(ip) = corr(curr_avg_fc,curr_spikes,'rows','pairwise','type','spearman');
end

%% Plot it
if 1
    figure
    pp = plot(spike_fc_corr,'o','linewidth',2);
    hold on
    plot(xlim,[0 0],'k--','linewidth',2)
    ylim([-1 1])
    median_corr = nanmedian(spike_fc_corr);
    xlabel('Patient')
    xticklabels([])
    ylabel('Spike rate - connectivity correlation (\rho)')
    [~,p] = ttest(spike_fc_corr);
    legend(pp,sprintf('median \\rho = %1.2f, %s',median_corr,get_p_text(p)),'location','southeast','fontsize',15)
    set(gca,'fontsize',15)
    print(gcf,[plot_folder,'spike_fc_corr'],'-dpng')
end

end