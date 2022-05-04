function nout = spike_fc_correlation

%% Parameters
do_plot = 0;
do_r2 = 0;
do_resid = 1;
max_spikes = 1/3600;

locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
model_folder = [results_folder,'analysis/outcome/plots/'];
atlas_folder = [results_folder,'analysis/atlas/'];


%% Load out file and get roc stuff
out = load([out_folder,'main_out.mat']);
out = out.out;

%% Also load spike out file
spout = load([out_folder,'spikes_out.mat']);
spout = spout.spikes_out;
sp_names = spout.names;
all_chs_corr = spout.all_chs_corr;
sp_chs_corr = spout.sp_chs_corr;

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

%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);

%% Take r^2
if do_r2
    fc = cellfun(@(x) x.^2, fc,'uniformoutput',false);
end

%% Do individual spike analysis
all_sp_chs_corr = nan(npts,2); % first column is pre spike and 2nd column is spike
all_all_chs_corr = nan(npts,2);
for ip = 1:npts
    curr_all = all_chs_corr{ip};
    curr_sp = sp_chs_corr{ip};
    
    if isempty(curr_all), continue; end
    
    % take mean across spikes
    curr_all_mean = nanmean(curr_all,1);
    curr_sp_mean = nanmean(curr_sp,1);
    
    all_sp_chs_corr(ip,:) =curr_sp_mean;
    all_all_chs_corr(ip,:) = curr_all_mean;
end
nout.all_sp_chs_corr = all_sp_chs_corr;
nout.all_all_chs_corr = all_all_chs_corr;

%% Spatially normalize the FC matrix
%[resid,f] = fit_distance_model(locs,fc,soz,rate,max_spikes,plot_folder);
dens_out = load([model_folder,'dens_model.mat']);
dens_out = dens_out.out;
resid = dens_out.resid;

%% Get spike-fc correlations
spike_fc_corr = nan(npts,1);
spike_resid_corr = nan(npts,1);
for ip = 1:npts
    curr_labels = labels{ip};
    assert(sum(find_non_intracranial(curr_labels)) == 0)
    
    curr_avg_fc = nanmean(fc{ip},2);
    curr_spikes = rate{ip};
    spike_fc_corr(ip) = corr(curr_avg_fc,curr_spikes,'rows','pairwise','type','spearman');
    
    curr_avg_resid = nanmean(resid{ip},2);
    spike_resid_corr(ip) = corr(curr_avg_resid,curr_spikes,'rows','pairwise','type','spearman');
end

nout.spike_fc_corr = spike_fc_corr;
nout.spike_resid_corr = spike_resid_corr;

%% Plot it
if do_plot
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