function univariate_fdr_plots

%% Parameters
which_pts = 'all';
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
fprintf(fid,'<p><br><b>Interictal EEG feature asymmetry better identifies left-sided SOZ than right-sided SOZ</b></br>');
fprintf(fid,['We first examined, using univariate comparisons, how interictal EEG feature AI differs '...
    'between patients with left-sided SOZs, right-sided SOZs, and bilateral SOZs.']);

%% Now do lr_mt to get AI features
[T,features] =  lr_mt(3); % just sleep 

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
set(gcf,'position',[1 1 1000 1000])
t = tiledlayout(2,2,"TileSpacing",'compact','padding','tight');


%% Show spikes
feature = 'spikes bipolar sleep';
%feature = 'plv theta bipolar wake';
nexttile
[~,stats] = boxplot_with_points(T.(feature),T.(response),1,{'left','right','bilateral'});
thresh = 0.4;
if 1
    plot(xlim,[thresh,thresh],'--','color',[0, 0.4470, 0.7410])
end
spikes = T.(feature);
left_soz_spikes = spikes(strcmp(T.(response),'left'));
right_soz_spikes = spikes(strcmp(T.(response),'right'));
ylabel({'Spike rate asymmetry index','(bipolar reference)'})
title('Spike rate asymmetry index by SOZ laterality')
set(gca,'fontsize',15)
%
fprintf(fid,[' Fig. 3A shows the mean spike rate AI (bipolar reference) '...
    'for patients with different SOZ lateralities. This feature was chosen '...
    'for visual example as the feature with the highest effect size at distinguish '...
    'SOZ lateralities. The spike rate AI differed across lateralities '...
    '(Kruskal-Wallis: &#967;<sup>2</sup><sub>%d</sub> = %1.1f, %s, '...
    '&#951;<sup>2</sup> = %1.2f). In post-hoc tests, the spike rate AI was higher in patients with '...
    'left-sided SOZs than in patients with right-sided SOZs (%s) and '...
    'in patients with bilateral SOZs (%s), but there was no difference '...
    'between patients with right and bilateral SOZs (%s). Patients with left-sided '...
    'SOZ tended to have positive spike rate AI, and patients with right-sided '...
    'SOZ tended to have negative spike rate AI, implying frequent concordance '...
    'between spike and seizure laterality. However, there was a high degree of overlap, '...
    'particularly between the right SOZ and bilateral SOZ group.</p>'],...
    stats.tbl{2,3},stats.tbl{2,5},get_p_html(stats.p),stats.eta2,...
    get_p_html(stats.lrp),get_p_html(stats.lbp),get_p_html(stats.rbp));
%}

if 1
    %% Example threshold - confusion matrix
    % Not cross validated, just see how the bipolar spike threshold (chosen
    % post hoc to be the best) separates left vs right/bilateral
    above_thresh = T.(feature) >thresh;
    below_thresh = T.(feature) < thresh;
    C(1,1) = sum(above_thresh & strcmp(T.(response),'left'));
    C(1,2) = sum(below_thresh & strcmp(T.(response),'left'));
    C(2,1) = sum(above_thresh & (strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')));
    C(2,2) = sum(below_thresh & (strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')));

    classes = {'Left','Right/bilateral'}; nclasses = 2;
    % Calculate accuracy
    accuracy = sum(diag(C))/sum(C(:));
    % Balanced accuracy is the average across all classes of the number of 
    % data accurately predicted belonging to class m divided by the number of
    % data belonging to class m
    recall = nan(nclasses,1);
    for i = 1:nclasses
        tp = C(i,i);
        fn = sum(C(i,~ismember(1:nclasses,i))); 
        recall(i) = tp/(tp+fn); % tp is number correctly predicted to be in class, tp + fn is everyone in the class
    end
    balanced_accuracy = mean(recall);
    sens_left = C(1,1)/(C(1,1)+C(1,2));
    spec_left = C(2,2)/(C(2,1)+C(2,2));
    

    nexttile(t)
    new_numbers = map_numbers_onto_range(C,[1 0]);
    Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
    D = diag(new_numbers);
    Dcolor = [repmat(D,1,2),ones(length(D),1)];
    Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
    imagesc(Ccolor)
    pretty_name = classes;
    xticks(1:nclasses)
    xticklabels((pretty_name))
    yticks(1:nclasses)
    yticklabels((pretty_name))
    ytickangle(90)
    xlabel('Predicted')
    ylabel('True')
    hold on
    for i = 1:nclasses
        for j = 1:nclasses
            text(i,j,sprintf('%d',C(j,i)),'horizontalalignment','center','fontsize',20)
        end
    end
    title(sprintf('Accuracy:%1.1f%%\nBalanced accuracy: %1.1f%%',...
        accuracy*100,balanced_accuracy*100))
    set(gca,'fontsize',15)

    fprintf(fid,[' As an example of how clinicians could immmediately use the AI'...
        ' to lateralize epilepsy, we built a simple decision rule using the spike'...
        ' rate AI (bipolar reference). Based on post-hoc analysis of Fig. 3A,'...
        ' we selected the AI that best appeared to separate left-sided SOZ from other '...
        'SOZ lateralities (0.4). We predicted that patients with a spike rate AI above 0.4 '...
        'would have left-sides SOZ, and patients with a spike rate AI below 0.4 would have right or '...
        'bilateral SOZ. Fig. 3B shows the resulting confusion matrix. The overall accuracy is %1.1f%%,'...
        ' the sensitivity to predict left-sided SOZ is %1.1f%%, and the specificity is %1.1f%%. '...
        'It is important to emphasize that this decision rule was created post-hoc and not cross-validated, '...
        'and thus should be interpreted with caution. However, it suggests that simply measuring average spike '...
        'rates in the left and right temporal lobes and calculating the asymmetry index may provide '...
        'helpful ancillary data to lateralize epilepsy.</p>'],accuracy*100,sens_left*100,spec_left*100);

else
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
        'most of the variance in the features (the number two was chosen for visualization purposes. '...
        'Fig. 3B shows the scores for the first two principal components for patients with '...
        'different SOZ lateralities. Patients with left-sided SOZs tend to cluster in the upper '...
        'right corner, with higher scores for both principal components. Patients with right-sided '...
        'SOZs tended to have lower scores for both principal components. Patients with bilateral SOZs '...
        'had scores centered around 0, overlapping with both unilateral groups.</p>']);
end

%% Univariate analyses of features
rm_sleep_text = @(x) strrep(x,' sleep','');
nfeatures = length(allowed_features);
feature_p_val = nan(nfeatures,1);
feature_eta2= nan(nfeatures,1);

% Loop over features
for i = 1:nfeatures

    % get p value and eta 2
    [feature_p_val(i),tbl] = kruskalwallis(T.(allowed_features{i}),T.(response),'off');
    feature_eta2(i) = tbl{2,2}/(tbl{2,2}+tbl{3,2});
    %feature_p_val(i) = anova1(T.(allowed_features{i}),T.(response),'off');
end
feature_p_val(isnan(feature_p_val)) = 1;

% False discovery rate
[~,qvalues] = mafdr(feature_p_val);
n_to_plot = 15;
[~,I] = sort(feature_eta2,'descend'); % sort by eta2
sorted_qvalues = qvalues(I(1:n_to_plot));

nexttile
plot(feature_eta2(I(1:n_to_plot)),'ko','markersize',15,'linewidth',2) % plot the top eta 2 scores
hold on
for i = 1:n_to_plot
    if sorted_qvalues(i) < 0.05
      %  plot(i,feature_eta2(I(i)),'k*','markersize',15,'linewidth',2)
    end
end
assert(all(sorted_qvalues<0.05)) % ensure q is actually < 0.05
%}
%ylim(yl_new)
xticks(1:n_to_plot)
ylabel('\eta^2_p')
xticklabels(cellfun(rm_sleep_text,cellfun(@greek_letters_plots,allowed_features(I(1:n_to_plot)),'uniformoutput',false),...
    'UniformOutput',false));
title('Effect size (\eta^2_p) to distinguish left/right/bilateral SOZ')
set(gca,'fontsize',15)

fprintf(fid,['<p>We next examined the ability of all interictal EEG features to '...
    'distinguish SOZ laterality. For each interictal EEG feature, we calculated the effect '...
    'size (&#951;<sup>2</sup>) at separating the three SOZ lateralities. We ranked '...
    'features in descending order by effect size. Fig. 3C shows the effect size of the top %d ranked features. '...
    'Each of these features had a significant effect at separating the three SOZ lateralities (Kruskal-Wallis test with '...
    'false discovery rate correction). Notably, many of the top-ranked AI '...
    'features involve spike rates, along with features related to '...
    'relative entropy, and bandpower. '],n_to_plot);

%% Univariate analyses for different comparisons
feature_p_val = nan(nfeatures,2);
feature_eta2= nan(nfeatures,2);
for i = 1:nfeatures
    curr_feat = T.(allowed_features{i});

    % Just separate each from bilateral
    if just_sep_bilat
        feature_p_val(i,1) = ranksum(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'bilateral')));
        feature_p_val(i,2) = ranksum(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'bilateral')));
    
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,1) = dT.Effect;
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,2) = dT.Effect;
    else
        feature_p_val(i,1) = ranksum(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')));
        feature_p_val(i,2) = ranksum(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')));
    
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'left')),curr_feat(strcmp(T.(response),'right')|strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,1) = dT.Effect;
        dT = meanEffectSize(curr_feat(strcmp(T.(response),'right')),curr_feat(strcmp(T.(response),'left')|strcmp(T.(response),'bilateral')),Effect="cohen");
        feature_eta2(i,2) = dT.Effect;
    end
end
feature_p_val(isnan(feature_p_val)) = 1;
% False discovery rate
[~,qvalues1] = mafdr(feature_p_val(:,1)); [~,qvalues2] = mafdr(feature_p_val(:,2));
n_to_plot = 15;
[~,I1] = sort(abs(feature_eta2(:,1)),'descend'); [~,I2] = sort(abs(feature_eta2(:,2)),'descend'); 
sorted_q1 = qvalues1(I1(1:n_to_plot)); sorted_q2 = qvalues2(I2(1:n_to_plot));
tt = tiledlayout(t,1,1);
tt.Layout.Tile = 4;
tt.Layout.TileSpan = [1 1];
ax1 = axes(tt);
ax2 = axes(tt);
ax2.XAxisLocation = 'bottom';
ax2.YAxisLocation = 'right';
pl = plot(ax1,abs(feature_eta2(I1(1:n_to_plot),1)),'o','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
hold on
pr = plot(ax2,abs(feature_eta2(I2(1:n_to_plot),2)),'o','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
%{
for i =1:n_to_plot
    if sorted_q1 < 0.05
        plot(ax1,i,abs(feature_eta2(I1(i),1)),'*','markersize',15,'color',[0 0.4470 0.7410],'linewidth',2);
    end

    if sorted_q2 < 0.05
        plot(ax2,i,abs(feature_eta2(I2(i),2)),'*','markersize',15,'color',[0.8500 0.3250 0.0980],'linewidth',2);
    end
end
%}
num_sig_1 = sum(sorted_q1<0.05);
num_sig_2 = sum(sorted_q2<0.05);
assert(num_sig_1 == n_to_plot && num_sig_2==0)
ax2.XColor = [0.8500 0.3250 0.0980];
ax2.YColor = [0.8500 0.3250 0.0980];
ax1.XColor = [0 0.4470 0.7410];
ax1.YColor = [0 0.4470 0.7410];
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
ax1.YLim = [min([min(abs(feature_eta2(I1(1:n_to_plot),1))),min(abs(feature_eta2(I2(1:n_to_plot),2)))])-0.2,...
    max([max(abs(feature_eta2(I1(1:n_to_plot),1))),max(abs(feature_eta2(I2(1:n_to_plot),2)))])+0.2];
ax2.YLim = [min([min(abs(feature_eta2(I1(1:n_to_plot),1))),min(abs(feature_eta2(I2(1:n_to_plot),2)))])-0.2,...
    max([max(abs(feature_eta2(I1(1:n_to_plot),1))),max(abs(feature_eta2(I2(1:n_to_plot),2)))])+0.2];
if 0
    table(feature_eta2(I2,2),qvalues2(I2),feature_p_val(I2,2))
end
ax1.XTick = 1:n_to_plot;
ax2.XTick = 1:n_to_plot;
ax1.XTickLabel = cellfun(rm_sleep_text,cellfun(@greek_letters_plots,allowed_features(I1(1:n_to_plot)),'uniformoutput',false),...
    'uniformoutput',false);
ax2.XTickLabel = cellfun(rm_sleep_text,cellfun(@greek_letters_plots,allowed_features(I2(1:n_to_plot)),'uniformoutput',false),...
    'uniformoutput',false);
ax1.XAxisLocation = 'top';


%xticklabels([])
if just_sep_bilat
    legend([pl pr],{'Left vs bilateral','Right vs bilateral'},'location','northeast','fontsize',15)
else
    legend([pl pr],{'Left vs right/bilateral','Right vs left/bilateral'},'location','northeast','fontsize',15)
end
%title('Effect sizes (Cohen''s {\it d}) to distinguish specific laterality')
set(ax1,'fontsize',15); set(ax2,'fontsize',15)
ylabel(ax1,'|Cohen''s {\it d}|','color','k','fontsize',15)

fprintf(fid,['We compared the ability of interictal feature AI to distinguish left-sided SOZs '...
    'versus right-sided SOZs. We separately calculated the absolute value of the effect size (Cohen''s <i>d</i>) '...
    'at distinguishing left-sided SOZs bilateral SOZs and that for distinguishing '...
    'right-sided SOZs from bilateral SOZs. For each of the two classification questions, '...
    'we ranked features in descending order by their absolute Cohen''s <i>d</i>. Fig. 3D '...
    'shows the Cohen''s <i>d</i> values for the top %d ranked features. '...
    'The top-ranked features are different for the two classification questions. For distinguishing '...
    'left-sided SOZ, spikes and relative entropy features performed best. For distinguishing '...
    'right-sided SOZ, several other features performed best. All of the top-%d ranked '...
    'features for distinguishing left-sided SOZs from right or bilateral SOZs had '...
    'significant effect sizes (Mann-Whitney U test with false discovery rate correction), '...
    'and <i>none</i> of the top %d ranked features at distinguishing right-sided SOZs from '...
    'left or bilateral SOZs had significant effect sizes after correcting for false '...
    'discovery rate.</p>'],n_to_plot);




print(gcf,[plot_folder,'Fig3'],'-dpng')

end