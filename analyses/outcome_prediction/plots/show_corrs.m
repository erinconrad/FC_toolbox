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
out = load([data_folder,'corr_out.mat']);
out = out.out;
avg_corr_sp = out.avg_corr_sp;
avg_corr_pear = out.avg_corr_pear;

%% Also load the atlas file
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;
atlas = out.atlas;
spikes = out.spikes_atlas;

%% Get normalized atlas
z = (atlas-nanmean(atlas,3))./nanstd(atlas,[],3);
ns = squeeze(nanmean(z,2));

%% Correlate ns and spikes
npts = size(atlas,3);
avg_corr_sp_norm = nan(npts,1);
for ip = 1:npts
    avg_corr_sp_norm(ip) = corr(ns(:,ip),spikes(:,ip),'rows','pairwise',...
        'type','spearman');
    
end

%% Initialize figure
figure
set(gcf,'position',[10 10 900 350])
tiledlayout(1,2,'tilespacing','tight','padding','tight')

%% Show spearman for non-normalized
nexttile
ind_corr_plot(avg_corr_sp)
title('Spike-connectivity correlation')

%% Show spearman for normalized
nexttile
ind_corr_plot(avg_corr_sp_norm)
title('Spike-connectivity correlation (normalized by anatomy)')


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