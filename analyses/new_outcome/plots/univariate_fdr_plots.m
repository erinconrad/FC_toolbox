function univariate_fdr_plots(T,features)

%% Parameters
response = 'soz_lats';

%{
which_sleep_stage = 'all';
allowed_features = features(cellfun(@(x) contains(x,which_sleep_stage),features));
%}
allowed_features = features;


%% Initialize figure
figure


%% Show spikes
feature = 'spikes delta car all';
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
ylim([0 0.6])
legend([pl pr],{'Left vs right/bilateral','Right vs left/bilateral'},'location','east','fontsize',15)
title('FDR-adjusted p-values to distinguish specific laterality')
set(gca,'fontsize',15)

%% Correlation matrix
%C = corr()

end