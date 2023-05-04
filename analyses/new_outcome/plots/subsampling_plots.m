function subsampling_plots

%% Parameters
pca_spikes_perc = 95;
rm_non_temporal = 1;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
if ~exist(plot_folder,'dir')
    mkdir(plot_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Run the lr_mt to extract features
[T,features,way,dur,sample,ss] =  lr_mt_multitime([2 3]);
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];

%% Establish things I am varying
all_ways = unique(way);
nways = length(all_ways);

all_durs = unique(dur);
ndurs = length(all_durs);

all_samples = unique(sample);
nsamples = length(all_samples);

all_ss = unique(ss);
nss = length(all_ss);


%% Prep figure
figure
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')

% Plan to save all data
all_data = nan(nss,2,ndurs,nsamples); % left, right;

% Loop over rand and continuous
for iw = 1
    curr_way = all_ways(iw);
    if curr_way == 1
        way_text = 'random';
    else
        way_text = 'continuous';
    end


    % Loop over nss
    for iss = 1:nss
        nexttile
        curr_ss = all_ss(iss);
        if curr_ss == 2
            ss_text = 'wake';
        elseif curr_ss == 3
            ss_text = 'sleep';
        end

        all_auc_l = nan(ndurs,nsamples);
        all_auc_r = nan(ndurs,nsamples);
        
        % Loop over ndurs -  these will be different error bar points
        for id = 1:ndurs
            fprintf('\nDoing way %d, ss %d, dur %d...',iw,iss,id);
            curr_dur = all_durs(id);
            
            
            % Loop over samples
            for is = 1:nsamples
                curr_sample = all_samples(is);

                % Get the relevant features
                relevant_features = contains(features,sprintf('_way%d_',curr_way)) & ...
                    contains(features,sprintf('_ss%d_',curr_ss)) & ...
                    contains(features,sprintf('_dur%d_',curr_dur)) & ...
                    contains(features,sprintf('_samp%d_',curr_sample));
                curr_features = features(relevant_features);

                % run the model
                just_spikes = 1; % Just spikes
                lefts = classifier_wrapper(T,curr_features,pca_spikes_perc,1,just_spikes,rm_non_temporal,[]);
                rights = classifier_wrapper(T,curr_features,pca_spikes_perc,2,just_spikes,rm_non_temporal,[]);

                % Get ROC stats
                [~,~,~,AUCL] = perfcurve(lefts.class,lefts.scores,lefts.pos_class);
                [~,~,~,AUCR] = perfcurve(rights.class,rights.scores,rights.pos_class);

                all_auc_l(id,is) = AUCL;
                all_auc_r(id,is) = AUCR;

            end
            fprintf('median left AUC is %1.2f\n',nanmedian(all_auc_l(id,:)))

        end

        % Get stats
        median_l = nanmedian(all_auc_l,2);
        median_r = nanmedian(all_auc_r,2);
        P_l_25 = prctile(all_auc_l,[25],2);
        P_r_25 = prctile(all_auc_r,[25],2);
        P_l_75 = prctile(all_auc_l,[75],2);
        P_r_75 = prctile(all_auc_r,[75],2);

        U_l = P_l_75-median_l;
        U_r = P_r_75-median_r;
        L_l = median_l - P_l_25;
        L_r = median_r - P_r_25;

        % Plot it
        el = errorbar(1:ndurs,median_l,L_l,U_l,'o','color',[0, 0.4470, 0.7410]);
        hold on
        er = errorbar(1:ndurs,median_r,L_r,U_r,'o','color',[0.8500, 0.3250, 0.0980]);
        ylim([0.4 1])

        legend([el,er],{'Left vs right/bilateral','Right vs left/bilateral'},'location','southeast')
        xticks(1:ndurs)
        xticklabels({'1 min','2 mins','5 mins','10 mins','20 mins'})
        ylabel('Median (IQR) AUC')
        title(sprintf('Model accuracy by duration\n(%s sampling, %s)',way_text,ss_text))
        set(gca,'fontsize',15)

        % save data
        all_data(iss,1,:,:) = all_auc_l;
        all_data(iss,2,:,:) = all_auc_r;

    end

end

print(gcf,[plot_folder,'subsample'],'-dpng')
savefig(gcf,[plot_folder,'subsample','.fig'])
save([plot_folder,'subsample.mat'],'all_data');

end