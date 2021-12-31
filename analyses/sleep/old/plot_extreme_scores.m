function plot_extreme_scores

%% parameters
ncomps = 2;
npts = 5;
which = 'seizure';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
int_folder = [results_folder,'analysis/intermediate/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


out = load([out_folder,'out.mat']);
out = out.out;
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

switch which
    case 'sleep'
        times = sleep_hist_out.times;
        bins = sleep_hist_out.all_pts_spikes_bins;
        names = sleep_hist_out.names;
    case 'seizure'
        times = sz_out.times;
        bins = sz_out.all_pts_spikes_bins;
        names = sz_out.names;
end

%% Normalize
norm_bins = (bins - nanmean(bins,2))./nanstd(bins,[],2);

%% PCA
[coeff,score,latent] = pca(norm_bins,'Rows','complete');

%% Prep figure
for ic = 1:ncomps
    figure
    set(gcf,'position',[10 10 1500 800])
    tiledlayout(2,5,'tilespacing','tight','padding','compact')

    for in = 1:2


        % sort pts by scores
        if in == 1 %top scores
            [sorted_score,I] = sort(score(:,ic),'descend');
            top_text = 'top';
        else
            [sorted_score,I] = sort(score(:,ic),'ascend');
            top_text = 'bottom';
        end
        I(isnan(sorted_score)) = [];

        for ip = 1:npts
            nexttile
            currp = I(ip);
            plot(times,bins(currp,:),'linewidth',2)
            hold on
            plot([0 0],ylim,'k--','linewidth',2)
            xlabel('Hours')
            ylabel('Spikes/elecs/min')
            title(sprintf('Component %d %s %d (%s)',ic,top_text,ip,names{I(ip)}))
            set(gca,'fontsize',15)
        end
    end

    print([out_folder,which,'_component_',num2str(ic),'_topscorers'],'-dpng')
    close(gcf)
end
    

end