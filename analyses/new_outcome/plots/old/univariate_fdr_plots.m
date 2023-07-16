function univariate_fdr_plots

%% Parameters
which_pts = 'hup';
rm_non_temporal = 1;
response = 'soz_lats';
just_sep_bilat = 1;

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
fprintf(fid,'<p><br><b>Several interictal EEG feature asymmetries are distinct across SOZ lateralities</b></br>');
fprintf(fid,['We compared interictal EEG features '...
    'between patients with left-sided, right-sided, and bilateral SOZs.']);

%% Now do lr_mt to get AI features
[T,features] =  lr_mt(3); % just sleep 
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

if 0
    feature = 'spikes bipolar sleep';
    figure
    boxplot_with_points(T.(feature),T.(response),0,{'left','right','bilateral'},[],'para');
    ylabel({'Spike rate asymmetry index','(bipolar reference)'})
    set(gca,'fontsize',25)
    print(gcf,[plot_folder,'methods_ex_AI'],'-depsc')

end

%% Initialize figure
figure
set(gcf,'position',[1 1 1000 1000])
t = tiledlayout(2,2,"TileSpacing",'compact','padding','tight');


%% A: Show spikes
feature = 'spikes bipolar sleep';
nexttile
[~,stats] = boxplot_with_points(T.(feature),T.(response),1,{'left','right','bilateral'},[],'para');
spikes = T.(feature);
ylabel({'Spike rate asymmetry index','(bipolar reference)'})
title('Spike rate asymmetry index by SOZ laterality')
set(gca,'fontsize',15)
%
fprintf(fid,[' Fig. 3A shows the mean spike rate AI (bipolar reference) '...
    'for patients with different SOZ lateralities. This feature was chosen '...
    'for visual example as the feature with the highest effect size at distinguish '...
    'SOZ lateralities. The spike rate AI differed across lateralities '...
    '(ANOVA: F(%d,%d) = %1.1f, %s, '...
    '&#951;<sup>2</sup> = %1.2f). In post-hoc tests, the spike rate AI was higher in patients with '...
    'left-sided SOZs than in patients with right-sided SOZs (%s) and '...
    'in patients with bilateral SOZs (%s), but there was no difference '...
    'between patients with right and bilateral SOZs after correcting for multiple comparisons (%s). Patients with left-sided '...
    'SOZ tended to have positive spike rate AI, and patients with right-sided '...
    'SOZ tended to have negative spike rate AI, implying frequent concordance '...
    'between spike and seizure laterality. However, there was a high degree of overlap, '...
    'particularly between the right SOZ and bilateral SOZ groups. '],...
    stats.tbl{2,3},stats.tbl{3,3},stats.tbl{2,5},get_p_html(stats.p),stats.eta2,...
    get_p_html(stats.lrp),get_p_html(stats.lbp),get_p_html(stats.rbp));
%}

%{
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
    nexttile(t)
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
        'most of the variance in the features (the number two was chosen for visualization purposes). '...
        'Fig. 3B shows the scores for the first two principal components for patients with '...
        'different SOZ lateralities. Patients with left-sided SOZs tend to cluster in the upper '...
        'right corner, with higher scores for both principal components. Patients with right-sided '...
        'SOZs tended to have lower scores for both principal components. Patients with bilateral SOZs '...
        'had scores centered around 0, overlapping with both unilateral groups.</p>']);
%}

%% Univariate analyses of features
rm_sleep_text = @(x) strrep(x,' sleep','');
nfeatures = length(allowed_features);
feature_p_val = nan(nfeatures,1);
feature_eta2= nan(nfeatures,1);

% Loop over features
for i = 1:nfeatures

    % get p value and eta 2
    [feature_p_val(i),tbl] = anova1(T.(allowed_features{i}),T.(response),'off');
    feature_eta2(i) = tbl{2,2}/(tbl{2,2}+tbl{3,2});
    
end
feature_p_val(isnan(feature_p_val)) = 1;

% False discovery rate
qvalues = mafdr(feature_p_val,'bhfdr',true); % Use BH approach, more conservative than Storey
n_to_plot = 15;
assert(sum(isnan(feature_eta2))==0)
[~,I] = sort(feature_eta2,'descend'); % sort by eta2
sorted_qvalues = qvalues(I(1:n_to_plot));

nexttile
plot(feature_eta2(I(1:n_to_plot)),'ko','markersize',15,'linewidth',2) % plot the top eta 2 scores
hold on
for i = 1:n_to_plot
    if sorted_qvalues(i) < 0.05
        plot(i,feature_eta2(I(i)),'k*','markersize',15,'linewidth',2)
    end
end
assert(all(sorted_qvalues<0.05)) % ensure q is actually < 0.05

xticks(1:n_to_plot)
ylabel('\eta^2_p')
xticklabels(cellfun(rm_sleep_text,cellfun(@greek_letters_plots,allowed_features(I(1:n_to_plot)),'uniformoutput',false),...
    'UniformOutput',false));
title('Effect size (\eta^2_p) to distinguish left/right/bilateral SOZ')
set(gca,'fontsize',15)

fprintf(fid,['We next examined the ability of all interictal EEG features to '...
    'distinguish SOZ laterality. For each interictal EEG feature, we calculated the effect '...
    'size (&#951;<sup>2</sup>) at separating the three SOZ lateralities. We ranked '...
    'features in descending order by effect size. Fig. 3B shows the effect size of the top %d ranked features. '...
    'Each of these features had a significant effect at separating the three SOZ lateralities (ANOVA with '...
    'Benjamini-Hochberg false discovery rate correction). The top-ranked AI '...
    'features involve spike rates, '...
    'relative entropy, and bandpower.</p>'],n_to_plot);

%% Univariate analyses for different comparisons
feature_p_val = nan(nfeatures,2);
feature_eta2= nan(nfeatures,2);
for i = 1:nfeatures
    curr_feat = T.(allowed_features{i});

    % Just separate each from bilateral
    if just_sep_bilat
        [~,feature_p_val(i,1)] = ttest2(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'bilateral')));
        [~,feature_p_val(i,2)] = ttest2(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'bilateral')));
    
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,1) = dT.Effect;
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,2) = dT.Effect;
    else
        [~,feature_p_val(i,1)] = ttest2(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')));
        [~,feature_p_val(i,2)] = ttest2(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')));
    
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,1) = dT.Effect;
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,2) = dT.Effect;
    end
end
%feature_p_val(isnan(feature_p_val)) = 1;
assert(sum(isnan(feature_eta2),'all')==0)
% False discovery rate
qvalues1 = mafdr(feature_p_val(:,1),'bhfdr',true); qvalues2 = mafdr(feature_p_val(:,2),'bhfdr',true);
n_to_plot = 15;
assert(sum(isnan(feature_eta2),'all')==0)
[~,I1] = sort(abs(feature_eta2(:,1)),'descend'); [~,I2] = sort(abs(feature_eta2(:,2)),'descend'); 
sorted_q1 = qvalues1(I1(1:n_to_plot)); sorted_q2 = qvalues2(I2(1:n_to_plot));
tt = tiledlayout(t,1,1);
tt.Layout.Tile = 3;
tt.Layout.TileSpan = [1 1];
ax1 = axes(tt);

pl = plot(ax1,1:n_to_plot,abs(feature_eta2(I1(1:n_to_plot),1)),'o','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
hold on
for i =1:n_to_plot
    if sorted_q1(i) < 0.05
        plot(ax1,i,abs(feature_eta2(I1(i),1)),'*','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
        
    end
end

ax1.XAxisLocation = 'top';
ax1.YAxisLocation = 'left';
ax1.XColor = [0 0.4470 0.7410];
ax1.YColor = [0 0.4470 0.7410];
ax1.Box = 'off';
ax1.YLim = [min([min(abs(feature_eta2(I1(1:n_to_plot),1))),min(abs(feature_eta2(I2(1:n_to_plot),2)))])-0.2,...
    max([max(abs(feature_eta2(I1(1:n_to_plot),1))),max(abs(feature_eta2(I2(1:n_to_plot),2)))])+0.2];
ax1.XTick = 1:n_to_plot;
ax1.XTickLabel = cellfun(rm_sleep_text,cellfun(@greek_letters_plots,allowed_features(I1(1:n_to_plot)),'uniformoutput',false),...
    'uniformoutput',false);


ax2 = axes(tt);
ax2.XAxisLocation = 'bottom';
ax2.YAxisLocation = 'right';
pr = plot(ax2,abs(feature_eta2(I2(1:n_to_plot),2)),'o','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
hold on
for i =1:n_to_plot
    if sorted_q2(i) < 0.05
        plot(ax2,i,abs(feature_eta2(I2(i),2)),'*','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
    end
end
ax2.XColor = [0.8500 0.3250 0.0980];
ax2.YColor = [0.8500 0.3250 0.0980];

ax2.Color = 'none';
ax2.Box = 'off';

ax2.YLim = [min([min(abs(feature_eta2(I1(1:n_to_plot),1))),min(abs(feature_eta2(I2(1:n_to_plot),2)))])-0.2,...
    max([max(abs(feature_eta2(I1(1:n_to_plot),1))),max(abs(feature_eta2(I2(1:n_to_plot),2)))])+0.2];
if 0
    table(feature_eta2(I2,2),qvalues2(I2),feature_p_val(I2,2))
end

ax2.XTick = 1:n_to_plot;
ax2.XTickLabel = cellfun(rm_sleep_text,cellfun(@greek_letters_plots,allowed_features(I2(1:n_to_plot)),'uniformoutput',false),...
    'uniformoutput',false);



%xticklabels([])
if just_sep_bilat
    legend([pl pr],{'Left vs bilateral','Right vs bilateral'},'location','northeast','fontsize',15)
else
    legend([pl pr],{'Left vs right/bilateral','Right vs left/bilateral'},'location','northeast','fontsize',15)
end
%title('Effect sizes (Cohen''s {\it d}) to distinguish specific laterality')
set(ax1,'fontsize',15); set(ax2,'fontsize',15)
ylabel(ax1,'|Cohen''s {\it d}|','color','k','fontsize',15)

fprintf(fid,['<p>We compared the set of interictal features that best distinguished '...
    'left from bilateral SOZs versus right from bilateral SOZs. '...
    'We separately calculated the absolute value of the effect size (Cohen''s <i>d</i>) '...
    'at distinguishing left-sided SOZs from bilateral SOZs and that for distinguishing '...
    'right-sided SOZs from bilateral SOZs. For each of the two classification questions, '...
    'we ranked features in descending order by their absolute Cohen''s <i>d</i>. Fig. 3C '...
    'shows the Cohen''s <i>d</i> values for the top %d ranked features. '...
    'For distinguishing '...
    'left-sided SOZ, spikes and relative entropy features performed best '...
    '(and two spike features were the only two that significantly discriminated left from bilateral). For distinguishing '...
    'right-sided SOZ, several other features performed best (though none were significant'...
    ' correcting for the false discovery rate).'],n_to_plot);

%% Concordance at the top plot
[conc_real,P,h] = alt_bootstrap_lr(T,features,1e3);
nexttile(t)
nfeatures = length(features);
chance_conc = ((1:nfeatures)/nfeatures)';
shaded_error_bars_fc(1:nfeatures,conc_real,P,'k')
hold on
pc = plot(1:nfeatures,chance_conc,'k--','linewidth',2);
ylabel('Proportion in common')
title('Top feature concordance for left versus right')
set(gca,'fontsize',15)
xlabel('Number of top-ranked features considered')
legend(pc,'Chance concordance','fontsize',15,'location','southeast')

assert(sum(h==1)==0)

fprintf(fid,[' We asked if there was a significant difference in the features that '...
    'best discriminated left from bilateral SOZs and those that best discriminated '...
    'right from bilateral SOZs. For every i ranging from 1 to the total number '...
    'of features (<i>N</i> = 180), we calculated the concordance in the top-i ranked '...
    'features according to Cohen''s <i>d</i> values for the two separate laterality '...
    'discriminations. We defined the concordance in the two feature sets to be '...
    '<i>Concordance</i><sub>i</sub> = length(intersection(<i>Features</i><sub>Left</sub>(i)'...
    ', <i>Features</i><sub>Right</sub>(i)))/<i>i</i>, where <i>Features</i><sub>Left</sub>(i)'...
    ' is the set of the top-i ranked features for discriminating left from '...
    'bilateral SOZs, and <i>Features</i><sub>Right</sub>(i) is the set of '...
    'the top-i ranked features for discriminating right from bilateral SOZs.'...
    ' We obtained 95%% confidence intervals (CI) for the concordance at each value of i'...
    ' by bootstrapping (<i>N</i> = 1,000), in which we randomly resampled patients'...
    ' with replacement and recalculated the concordance. We compared the 95%% CI'...
    ' to the chance concordance, which equals <i>concordance</i><sub>chance</sub>(i) = '...
    ' <i>i/N</i>, where <i>N</i> is the total number of features. There was no value '...
    'of i for which the chance concordance fell outside the 95%% CI, implying that '...
    'the best features to identify left-sided SOZs are not clearly different from the best '...
    'features to identify right-sided SOZs (nor are they more similar than expected by chance) (Fig. 3D). </p>']);

print(gcf,[plot_folder,'Fig3'],'-dpng')




end