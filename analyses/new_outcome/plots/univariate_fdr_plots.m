function univariate_fdr_plots

%% Parameters
rm_non_temporal = 0;
response = 'soz_lats';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];
subplot_path = [plot_folder,'ai_subplots/'];
if ~exist(subplot_path,'dir')
    mkdir(subplot_path)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Now do lr_mt to get AI features
[T,features] =  lr_mt; 

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
title('Spike rate asymmetry index by SOZ laterality')
set(gca,'fontsize',15)

%% Univariate analyses of features
nfeatures = length(allowed_features);
feature_p_val = nan(nfeatures,1);
feature_eta2= nan(nfeatures,1);
for i = 1:nfeatures
    [feature_p_val(i),tbl] = kruskalwallis(T.(allowed_features{i}),T.(response),'off');
    feature_eta2(i) = tbl{2,2}/(tbl{2,2}+tbl{3,2});
    %feature_p_val(i) = anova1(T.(allowed_features{i}),T.(response),'off');
end
feature_p_val(isnan(feature_p_val)) = 1;

% False discovery rate
[~,qvalues] = mafdr(feature_p_val);
n_to_plot = 20;
[~,I] = sort(feature_eta2,'descend');
sorted_qvalues = qvalues(I(1:n_to_plot));

nexttile
plot(feature_eta2(I(1:n_to_plot)),'ko','markersize',15,'linewidth',2)
hold on
for i = 1:n_to_plot
    if sorted_qvalues(i) < 0.05
        plot(i,feature_eta2(I(i)),'k*','markersize',15,'linewidth',2)
    end
end
%}
%ylim(yl_new)
xticks(1:n_to_plot)
ylabel('\eta^2_p')
xticklabels(cellfun(@greek_letters_plots,allowed_features(I(1:n_to_plot)),'uniformoutput',false))
title('Effect size (\eta^2_p) to distinguish left/right/bilateral SOZ')
set(gca,'fontsize',15)

%% Univariate analyses for different comparisons
feature_p_val = nan(nfeatures,2);
feature_eta2= nan(nfeatures,2);
for i = 1:nfeatures
    curr_feat = T.(allowed_features{i});
    feature_p_val(i,1) = ranksum(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')));
    feature_p_val(i,2) = ranksum(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')));

    dT = meanEffectSize(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')),Effect="cohen");
    feature_eta2(i,1) = dT.Effect;
    dT = meanEffectSize(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')),Effect="cohen");
    feature_eta2(i,2) = dT.Effect;
end
feature_p_val(isnan(feature_p_val)) = 1;
% False discovery rate
[~,qvalues1] = mafdr(feature_p_val(:,1)); [~,qvalues2] = mafdr(feature_p_val(:,2));
n_to_plot = 20;
[~,I1] = sort(abs(feature_eta2(:,1)),'descend'); [~,I2] = sort(abs(feature_eta2(:,2)),'descend'); 
sorted_q1 = qvalues1(I1(1:n_to_plot)); sorted_q2 = qvalues2(I2(1:n_to_plot));
nexttile
pl = plot(abs(feature_eta2(I1(1:n_to_plot),1)),'o','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
hold on
pr = plot(abs(feature_eta2(I2(1:n_to_plot),2)),'o','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
for i =1:n_to_plot
    if sorted_q1 < 0.05
        plot(i,abs(feature_eta2(I1(i),1)),'*','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
    end

    if sorted_q2 < 0.05
        plot(i,abs(feature_eta2(I2(i),2)),'*','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
    end
end

if 0
    table(feature_eta2(I2,2),qvalues2(I2),feature_p_val(I2,2))
end

xticklabels([])
legend([pl pr],{'Left vs right/bilateral','Right vs left/bilateral'},'location','northeast','fontsize',15)
title('Effect sizes (Cohen''s {\it d}) to distinguish specific laterality')
set(gca,'fontsize',15)
xlabel('Features (ordered by effect size)')
ylabel('|Cohen''s {\it d}|')


%% PCA/clustering 
% need to decide about imputation, etc.
X = table2array(T(:,features));
soz_lats = T.soz_lats;
% for now just remove rows with any nans
nan_rows = any(isnan(X),2);
X(nan_rows,:) = [];
soz_lats(nan_rows) = [];
right = strcmp(soz_lats,'right');
left = strcmp(soz_lats,'left');
bilateral = strcmp(soz_lats,'bilateral');

% normalize
X = (X-nanmean(X,1))./nanstd(X,[],1);
[coeff,score,latent] = pca(X);
nexttile
plot(score(left,1),score(left,2),'o','markersize',12,'linewidth',2)
hold on
plot(score(right,1),score(right,2),'+','markersize',12,'linewidth',2)
plot(score(bilateral,1),score(bilateral,2),'*','markersize',12,'linewidth',2)
xlabel('Component 1 score')
ylabel('Component 2 score')
title('Feature separation by SOZ laterality')
legend({'left','right','bilateral'},'location','southeast','fontsize',15)
set(gca,'fontsize',15)

print(gcf,[plot_folder,'Fig3'],'-dpng')

end