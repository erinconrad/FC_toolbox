function show_corrs

%% Parameters
which_atlas = 'brainnetome';%'aal_bernabei';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/plots/'];
data_folder = [results_folder,'analysis/outcome/data/'];
atlas_folder = [results_folder,'analysis/atlas/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the corr_out file
out = load([data_folder,'main_out.mat']);
out = out.out;
soz = out.all_soz_bin;
fc = out.all_fc;
locs = out.all_locs;
alt_spikes = out.all_spikes;

%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);

%% Spatially normalize the FC matrix
resid = fit_distance_model(locs,fc,soz);

%% Get ns of fc and resid
alt_ns = cellfun(@(x) nanmean(x,2),fc,'uniformoutput',false);
ns_resid = cellfun(@(x) nanmean(x,2),resid,'uniformoutput',false);

%% Also load the atlas file
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;
atlas = out.atlas;
spikes = out.spikes_atlas;

%% Get normalized atlas
z = (atlas-nanmean(atlas,3))./nanstd(atlas,[],3);
z = squeeze(nanmean(z,2));
ns = squeeze(nanmean(atlas,2));

%% Correlations
npts = size(atlas,3);
norm_ana = nan(npts,1);
base = nan(npts,1);
base_alt = nan(npts,1);
norm_dist = nan(npts,1);
for ip = 1:npts
    norm_ana(ip) = corr(z(:,ip),spikes(:,ip),'rows','pairwise',...
        'type','spearman');
    base(ip) = corr(ns(:,ip),spikes(:,ip),'rows','pairwise',...
        'type','spearman');
    base_alt(ip) = corr(alt_ns{ip},alt_spikes{ip},'rows','pairwise',...
        'type','spearman');
    norm_dist(ip) = corr(ns_resid{ip},alt_spikes{ip},'rows','pairwise',...
        'type','spearman');
    
end

%% Initialize figure
figure
set(gcf,'position',[10 10 900 350])
tiledlayout(1,5,'tilespacing','tight','padding','tight')

%% Show spearman for non-normalized
nexttile
ind_corr_plot(base)
title('Spike-connectivity correlation (atlas file)')

%% Show spearman for normalized
nexttile
ind_corr_plot(base_alt)
title('Spike-connectivity correlation (main file)')

nexttile
ind_corr_plot(norm_dist)
title('Spike-connectivity correlation (distance norm)')

nexttile
ind_corr_plot(norm_ana)
title('Spike-connectivity correlation (ana norm)')


end

function ind_corr_plot(corr_thing)

%% parameters
colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980];

plot(corr_thing,'o','linewidth',2,'color',colors(1,:))
hold on
plot(xlim,[nanmedian(corr_thing) nanmedian(corr_thing)],'linewidth',2,'color',colors(1,:))
plot(xlim,[0 0],'k--','linewidth',2)
ylim([-1 1])

% one sample ttest
[~,p,~,stats] = ttest(corr_thing);
xl = xlim;
yl = ylim;
text(xl(1),yl(2),sprintf('%s',get_p_text(p)),'verticalalignment','top','fontsize',15)
set(gca,'fontsize',15)
ylabel('Correlation coefficient')
xlabel('Patient')
xticklabels([])

end