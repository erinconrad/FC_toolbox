function univariate_fdr_plots(T,features)

%% Parameters
rm_non_temporal = 1;
response = 'soz_lats';

%{
which_sleep_stage = 'all';
allowed_features = features(cellfun(@(x) contains(x,which_sleep_stage),features));
%}
allowed_features = features;

%% Remove non temporal patients
if rm_non_temporal
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];
end

%% Initialize figure
figure
set(gcf,'position',[1 1 1400 1000])
tiledlayout(2,2,"TileSpacing",'tight','padding','tight')


%% Show spikes
feature = 'spikes car sleep';
%feature = 'plv theta bipolar wake';
nexttile
boxplot_with_points(T.(feature),T.(response),1)
ylabel('Spike rate asymmetry index')
set(gca,'fontsize',15)

%% Univariate analyses of features
nfeatures = length(allowed_features);
feature_p_val = nan(nfeatures,1);
for i = 1:nfeatures
    feature_p_val(i) = kruskalwallis(T.(allowed_features{i}),T.(response),'off');
end
feature_p_val(isnan(feature_p_val)) = 1;

% False discovery rate
[~,qvalues] = mafdr(feature_p_val);
n_to_plot = 20;
[~,I] = sort(qvalues,'ascend');

nexttile
plot(qvalues(I(1:n_to_plot)),'ko','markersize',15)
hold on
plot(xlim,[0.05 0.05],'k--')
xticks(1:n_to_plot)
xticklabels(cellfun(@greek_letters_plots,allowed_features(I(1:n_to_plot)),'uniformoutput',false))
title('FDR-adjusted p-values to distinguish left/right/bilateral SOZ')
ylim([0 0.07])
set(gca,'fontsize',15)

%% Univariate analyses for different comparisons
feature_p_val = nan(nfeatures,2);
for i = 1:nfeatures
    curr_feat = T.(allowed_features{i});
    feature_p_val(i,1) = ranksum(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')));
    feature_p_val(i,2) = ranksum(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')));
end
feature_p_val(isnan(feature_p_val)) = 1;
% False discovery rate
[~,qvalues1] = mafdr(feature_p_val(:,1)); [~,qvalues2] = mafdr(feature_p_val(:,2));
n_to_plot = 20;
[~,I1] = sort(qvalues1,'ascend'); [~,I2] = sort(qvalues2,'ascend'); 
nexttile
pl = plot(qvalues1(I1(1:n_to_plot)),'o','markersize',15);
hold on
pr = plot(qvalues2(I2(1:n_to_plot)),'o','markersize',15);
plot(xlim,[0.05 0.05],'k--')
xticklabels([])
ylim([0 0.3])
legend([pl pr],{'Left vs right/bilateral','Right vs left/bilateral'},'location','northeast','fontsize',15)
title('FDR-adjusted p-values to distinguish specific laterality')
set(gca,'fontsize',15)

%% For all/CAR, how do different AIs agree with eachother
sleep_stage = 'sleep';
ref = 'car';
%networks = {'pearson','coh','plv','re','xcor'};
curr_features = features(cellfun(@(x) contains(x,sleep_stage) && contains(x,ref) && ~contains(x,'SD'),features));
%matofindices = unique(cell2mat(arrayfun(@(pat)(find(contains(curr_features,pat))),networks(:),'UniformOutput',false)));
%curr_networks = curr_features(matofindices);
C = corr(table2array(T(:,curr_features)),'rows','pairwise');

nexttile
turn_nans_gray(C)
xticks(1:length(curr_features))
yticks(1:length(curr_features))
xticklabels(curr_features)
yticklabels(curr_features)
colorbar
clim([-1 1])
title('Correlation between AI features')
set(gca,'fontsize',15)

end