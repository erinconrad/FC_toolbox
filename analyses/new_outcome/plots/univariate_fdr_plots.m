function univariate_fdr_plots

%% Parameters
which_pts = 'hup';
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

%% Initialize results file
fname = [plot_folder,'results.html'];
fid = fopen(fname,'a');
fprintf(fid,'<p><br><b>Interictal EEG feature asymmetry better identifies left-sided SOZ than right-sided SOZ</b></br>');
fprintf(fid,['We first examined, using univariate comparisons, how interictal EEG feature AI differs '...
    'between patients with left-sided SOZs, right-sided SOZs, and bilateral SOZs.']);

%% Now do lr_mt to get AI features
[T,features] =  lr_mt; 

%{
which_sleep_stage = 'all';
allowed_features = features(cellfun(@(x) contains(x,which_sleep_stage),features));
%}
allowed_features = features;

%% Restrict to desired hospital
switch which_pts
    case 'all'
    case 'hup'
        hup = contains(T.names,'HUP');
        T(~hup,:) = [];
    case 'musc'
        musc = contains(T.names,'MP');
        T(~musc,:) = [];
end

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
[~,stats] = boxplot_with_points(T.(feature),T.(response),1,{'left','right','bilateral'});
spikes = T.(feature);
left_soz_spikes = spikes(strcmp(T.(response),'left'));
right_soz_spikes = spikes(strcmp(T.(response),'right'));
ylabel('Spike rate asymmetry index')
title('Spike rate asymmetry index by SOZ laterality')
set(gca,'fontsize',15)
%
fprintf(fid,[' Fig. 3A shows the spike rate AI (in sleep, common average reference) '...
    'for patients with different SOZ lateralities. The spike rate AI differed across lateralities '...
    '(Kruskal-Wallis: &#967;<sup>2</sup><sub>%d</sub> = %1.1f, %s, '...
    '&#951;<sup>2</sup> = %1.2f). In post-hoc tests, the spike rate AI was higher in patients with '...
    'left-sided SOZs than in patients with right-sided SOZs (%s) and '...
    'in patients with bilateral SOZs (%s), but there was no difference '...
    'between patients with right and bilateral SOZs (%s). In %d of %d patients '...
    'with left-sided SOZs, the spike rate AI was positive, as expected. However, in %d of %d '...
    'patients with right-sided SOZs, the spike rate AI was positive, implying frequent discordance '...
    'between temporal lobe spike rate laterality and ictal laterality for patients with right-sided SOZs.</p>'],...
    stats.tbl{2,3},stats.tbl{2,5},get_p_html(stats.p),stats.eta2,...
    get_p_html(stats.lrp),get_p_html(stats.lbp),get_p_html(stats.rbp),...
    sum(left_soz_spikes>0),length(left_soz_spikes),...
    sum(right_soz_spikes>0),length(right_soz_spikes));
%}
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

fprintf(fid,['<p>We next examined the ability of all interictal EEG features to '...
    'distinguish SOZ laterality. For each interictal EEG feature, we calculated the effect '...
    'size (&#951;<sup>2</sup>) at separating the three SOZ lateralities. We ranked '...
    'features in descending order by effect size. Fig. 3B shows the effect size of the top 20 ranked features. '...
    'Each of these features had a significant effect at separating the three SOZ lateralities (Kruskal-Wallis test with '...
    'false discovery rate correction). Notably, many of the top-ranked AI '...
    'features involve spike rates, along with features related to '...
    'relative entropy, Pearson correlation, bandpower, and coherence. ']);

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

fprintf(fid,['We compared the ability of interictal feature AI to distinguish left-sided SOZs '...
    'versus right-sided SOZs. We separately calculated the absolute value of the effect size (Cohen''s <i>d</i>) '...
    'at distinguishing left-sided SOZs from right or bilateral SOZs and that for distinguishing '...
    'right-sided SOZs from left or bilateral SOZs. For each of the two classification questions, '...
    'we ranked features in descending order by their absolute Cohen''s <i>d</i>. Fig. 3C '...
    'shows the Cohen''s <i>d</i> values for the top-20 ranked features, '...
    '(the top-ranked features are different for the two classification questions). All of the top-20 ranked '...
    'features for distinguishing left-sided SOZs from right or bilateral SOZs had '...
    'significant effect sizes (Mann-Whitney U test with false discovery rate correction), '...
    'and none of the top-20 ranked features at distinguishing right-sided SOZs from '...
    'left or bilateral SOZs had significant effect sizes after correcting for false '...
    'discovery rate.</p>']);


%% PCA/clustering 
% need to decide about imputation, etc.
X = table2array(T(:,features));
soz_lats = T.soz_lats;
%{
% for now just remove rows with any nans
nan_rows = any(isnan(X),2);
X(nan_rows,:) = [];
soz_lats(nan_rows) = [];
%}

% Imputation of nans. Make Nans equal to mean across other rows for that
% column
for ic = 1:size(X,2)
    nan_row = isnan(X(:,ic));
    X(nan_row,ic) = nanmean(X(:,ic));
end

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

fprintf(fid,['<p>We next attempted to visualize the ability of the set of interictal EEG '...
    'AI features as a whole to separate epilepsy lateralities given the high '...
    'degree of inter-feature correlation. We performed PCA across all features after '...
    'normalizing the features by subtracting the mean and dividing by the standard '...
    'deviation across patients. We obtained the two principal components that explained '...
    'most of the variance in the features (the number two was chosen for visualization purposes. '...
    'Fig. 3D shows the scores for the first two principal components for patients with '...
    'different SOZ lateralities. Patients with left-sided SOZs tend to cluster in the upper '...
    'right corner, with higher scores for both principal components. However, there is no '...
    'clear separation in the first two principal components between patients with right-sided '...
    'SOZs and patients with bilateral SOZs</p>.']);

print(gcf,[plot_folder,'Fig3'],'-dpng')

end