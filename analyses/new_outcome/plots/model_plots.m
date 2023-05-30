function model_plots


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

%% Load the intermediat file
out = load([plot_folder,'ext_models.mat']);
out = out.all;
model = out.model;
nmodels = length(model);

%% Initialize results file
fname = [plot_folder,'results.html'];
fid = fopen(fname,'a');
fprintf(fid,'<p><br><b>A classifier incorporating interictal EEG features predicts SOZ laterality</b></br>');


%% Initialize figure
figure
set(gcf,'position',[1 1 1200 1400])
tiledlayout(3,3,"TileSpacing",'tight','padding','tight')

%% A-C Do ROC curves for internal validation
for iv = 1
    for im = 1:nmodels
        nexttile
        ll = plot(model(im).val(iv).side(1).result.X,...
            model(im).val(iv).side(1).result.Y,'linewidth',2);
        hold on
        lr = plot(model(im).val(iv).side(2).result.X,...
            model(im).val(iv).side(2).result.Y,'linewidth',2);
        plot([0 1],[0 1],'k--','linewidth',2)
        xlabel('False positive rate')
        ylabel('True positive rate')
        legend([ll,lr],{sprintf('%s: AUC = %1.2f',...
            model(im).val(iv).side(1).description,model(im).val(iv).side(1).result.AUC),...
            sprintf('%s: AUC = %1.2f',...
            model(im).val(iv).side(2).description,model(im).val(iv).side(2).result.AUC)},'fontsize',15,...
            'location','southeast')
        title(sprintf('%s %s',model(im).type,model(im).val(iv).description))
        set(gca,'fontsize',15)


    end
end

fprintf(fid,['We tested whether interictal features can predict SOZ laterality in unseen patients. '...
    'We first performed leave-one-out classification on the HUP patients using two'...
    ' separate logistic regression models: one classifying a patient as having either a left-sided SOZ '...
    'as opposed to a bilateral or right-sided SOZ, and the other classifying right versus left or'...
    ' bilateral. We first tested these models using all interictal features, with PCA for dimensionality '...
    'reduction and LASSO for regularization and feature selection (see Methods). '...
    'The AUC of the ROC of a model predicting left vs right/bilateral SOZs using all features was %1.2f, '...
    'and one predicting right vs left/bilateral SOZs was %1.2f (Fig. 5A). We compared this to a model '...
    'trained only on spike rates (bipolar reference), whose corresponding AUCs were %1.2f '...
    'and %1.2f for the left and right models, respectively (Fig. 5B). Finally, a model trained only on '...
    'spike rates, where we converted each spike rate asymmetry index into a binary variable '...
    '(1=more spikes on the left, 0=more spikes on the right) achieved an AUC of %1.2f and %1.2f '...
    'for the left and right models, respectively (Fig. 5C). These results suggest that a model incorporating '...
    'many interictal features lateralizes temporal lobe epilepsy with moderate accuracy. '...
    'A model using only spike rates outperforms this complex model, although converting spike rates '...
    'to a simple binary variable indicating which side has more spikes dramatically reduces performance.</p>'],...
    model(1).val(1).side(1).result.AUC,model(1).val(1).side(2).result.AUC,...
    model(2).val(1).side(1).result.AUC,model(2).val(1).side(2).result.AUC,...
    model(3).val(1).side(1).result.AUC,model(3).val(1).side(2).result.AUC);

%% D: Subsampling plots, internal validation
if 1
sub = out.cv_ss; % cross validation

nexttile

durations = [1 5 10 20 30];
curr_ss = 2; % just do sleep
ndurs = length(durations);
nsamples = size(sub,4);


auc_l = squeeze(sub(curr_ss,1,:,:));
auc_r = squeeze(sub(curr_ss,2,:,:));

median_l = nanmedian(auc_l,2);
median_r = nanmedian(auc_r,2);
P_l_25 = prctile(auc_l,[25],2);
P_r_25 = prctile(auc_r,[25],2);
P_l_75 = prctile(auc_l,[75],2);
P_r_75 = prctile(auc_r,[75],2);

U_l = P_l_75-median_l;
U_r = P_r_75-median_r;
L_l = median_l - P_l_25;
L_r = median_r - P_r_25;

% Plot it


el = shaded_error_bars_fc(1:ndurs,median_l,[P_l_75';P_l_25'],[0, 0.4470, 0.7410]);
hold on
er = shaded_error_bars_fc(1:ndurs,median_r,[P_r_75';P_r_25'],[0.8500, 0.3250, 0.0980]);


errorbar(1:ndurs,median_l,L_l,U_l,'o','color',[0, 0.4470, 0.7410],...
    'LineWidth',2,'MarkerSize',10);
hold on
errorbar(1:ndurs,median_r,L_r,U_r,'o','color',[0.8500, 0.3250, 0.0980],...
    'LineWidth',2,'MarkerSize',10);
%}

ylim([0.55 0.9])

legend([el,er],{'Left vs right/bilateral','Right vs left/bilateral'},...
    'location','southeast','fontsize',15)
xticks(1:ndurs)
xticklabels(arrayfun(@(x) sprintf('%d min',x),durations,'uniformoutput',false))
ylabel('Median (IQR) AUC')
title(sprintf('Spike model accuracy by duration\n(HUP cross-validation)'))
set(gca,'fontsize',15)

fprintf(fid,['<p>Given the high performance of the spike rate-only model, and given '...
    'that this feature could easily be calculated in clinical settings, we performed '...
    'additional analyses to probe the accuracy of this model. We first examined how the '...
    'model accuracy depends on the length of time sampled. We selected random one-minute '...
    'segments amounting to several different total durations (1 minute, 5 minutes, '...
    '10 minutes, 20 minutes, and 30 minutes) and reconstructed the models, using ten different '...
    'random resamplings to generate statistics. Fig. 5D shows the median and IQR AUCs for both '...
    'the left and the right models. Model accuracies rise quickly with duration sampled, achieving '...
    'an accuracy similar to the full-duration models with 10 minutes of sampling.']);

end

%% E-F Confusion matrix (spikes only, internal cross-validation)
curr = model(2).val(1);
both_bal_acc = nan(2,1);
for is = 1:2
    scores = curr.side(is).result.scores;
    class = curr.side(is).result.class;

    classes = curr.side(is).result.unique_classes;
    nclasses = length(classes);
    pos_class = curr.side(is).result.pos_class;
    neg_class = classes(~strcmp(classes,pos_class));

    % put the positive class on top
    if ~strcmp(curr.side(is).result.pos_class,classes{1})
        classes = classes([2 1]);
    end

    % Get optimal ROC point
    [X,Y,T,~,opt] = perfcurve(class,scores,pos_class);
    %opt_thresh = T((X==opt(1))&(Y==opt(2)));
    opt_thresh = my_opt(X,Y,T);

    
    pred = cell(length(class),1);
    pred(scores >= opt_thresh) = {pos_class};
    pred(scores < opt_thresh) = neg_class;

    % Rebuild C
    positive = strcmp(class,pos_class);
    negative = strcmp(class,neg_class);
    pred_positive = strcmp(pred,pos_class);
    pred_negative = strcmp(pred,neg_class);
    C(1,1) = sum(positive & pred_positive);
    C(1,2) = sum(positive & pred_negative);
    C(2,1) = sum(negative & pred_positive);
    C(2,2) = sum(negative & pred_negative);
    
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
    both_bal_acc(is) = balanced_accuracy;
    
    % Plot
    nexttile
    % Map numbers onto 0 to 1
    new_numbers = map_numbers_onto_range(C,[1 0]);
    Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
    D = diag(new_numbers);
    Dcolor = [repmat(D,1,2),ones(length(D),1)];
    Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
    imagesc(Ccolor)
    
    % replace classnames
    pretty_name = classes;
    pretty_name = strrep(pretty_name,'left','Left');
    pretty_name = strrep(pretty_name,'right','Right');
    pretty_name = strrep(pretty_name,'br','Right/bilateral');
    pretty_name = strrep(pretty_name,'bl','Left/bilateral');
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
    title(sprintf('Spikes (HUP cross-validation)\nBalanced accuracy %1.1f%%',balanced_accuracy*100))
    set(gca,'fontsize',15)

end

fprintf(fid,[' To demonstrate statistics for a specific use-case, we selected the '...
    'operating points of the spike-model ROC curves with the highest summed true positive '...
    'rate and (1-false positive rate). Confusion matrices for the left- and right-sided models at these '...
    'operating points are shown in Figs. 5E and 5F. The balanced accuracy '...
    'was %1.1f%% for the model predicting left vs. right/bilateral SOZ, '...
    'and %1.1f%% for the model predicting right vs. left/bilateral SOZ.</p>'],...
    both_bal_acc(1)*100,both_bal_acc(2)*100);

%% G: ROC for spikes, external validation
iv = 2;
im = 2;
nexttile
ll = plot(model(im).val(iv).side(1).result.X,...
    model(im).val(iv).side(1).result.Y,'linewidth',2);
hold on
lr = plot(model(im).val(iv).side(2).result.X,...
    model(im).val(iv).side(2).result.Y,'linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr],{sprintf('%s: AUC = %1.2f',...
    model(im).val(iv).side(1).description,model(im).val(iv).side(1).result.AUC),...
    sprintf('%s: AUC = %1.2f',...
    model(im).val(iv).side(2).description,model(im).val(iv).side(2).result.AUC)},...
    'location','southeast','fontsize',15)
title(sprintf('%s %s',model(im).type,model(im).val(iv).description))
set(gca,'fontsize',15)

fprintf(fid,['<p>Finally, we tested how the spike-only models performed in an '...
    'external dataset. We trained the models using only spike rate data '...
    'on all HUP patients, and tested the models on patients from MUSC. '...
    'The model AUCs on the MUSC patients were %1.2f and %1.2f for '...
    'the left-sided and right-sided models, respectively (Fig. 5G).'],...
    model(2).val(2).side(1).result.AUC,model(2).val(2).side(2).result.AUC);


%% H-I: Confusion matrix (spikes only, external validation)
curr = model(2).val(2); % spike  model, external validation
both_bal_acc = nan(2,1);
for is = 1:2
    scores = curr.side(is).result.scores;
    class = curr.side(is).result.class;

    classes = curr.side(is).result.unique_classes;
    nclasses = length(classes);
    pos_class = curr.side(is).result.pos_class;
    neg_class = classes(~strcmp(classes,pos_class));

    % put the positive class on top
    if ~strcmp(curr.side(is).result.pos_class,classes{1})
        classes = classes([2 1]);
    end

    % Get optimal ROC point
    [X,Y,T,~,opt] = perfcurve(class,scores,pos_class);
    %opt_thresh = T((X==opt(1))&(Y==opt(2)));
    opt_thresh = my_opt(X,Y,T);
    
    pred = cell(length(class),1);
    pred(scores >= opt_thresh) = {pos_class};
    pred(scores < opt_thresh) = neg_class;

    % Rebuild C
    positive = strcmp(class,pos_class);
    negative = strcmp(class,neg_class);
    pred_positive = strcmp(pred,pos_class);
    pred_negative = strcmp(pred,neg_class);
    C(1,1) = sum(positive & pred_positive);
    C(1,2) = sum(positive & pred_negative);
    C(2,1) = sum(negative & pred_positive);
    C(2,2) = sum(negative & pred_negative);
    
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
    both_bal_acc(is) = balanced_accuracy;
    
    % Plot
    nexttile
    % Map numbers onto 0 to 1
    new_numbers = map_numbers_onto_range(C,[1 0]);
    Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
    D = diag(new_numbers);
    Dcolor = [repmat(D,1,2),ones(length(D),1)];
    Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
    imagesc(Ccolor)
    
    % replace classnames
    pretty_name = classes;
    pretty_name = strrep(pretty_name,'left','Left');
    pretty_name = strrep(pretty_name,'right','Right');
    pretty_name = strrep(pretty_name,'br','Right/bilateral');
    pretty_name = strrep(pretty_name,'bl','Left/bilateral');
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
    title(sprintf('Spikes (MUSC external validation)\nBalanced accuracy %1.1f%%',balanced_accuracy*100))
    set(gca,'fontsize',15)

end

fprintf(fid,[' The corresponding confusion matrices for the MUSC patients are '...
    'shown in Fig. 5H and 5I, with balanced accuracies of %1.1f%% and %1.1f%% '...
    'for the left-sided and right-sided models, respectively.'],both_bal_acc(1)*100,...
    both_bal_acc(2)*100);

fprintf(fid,[' These results suggest that models '...
    'using only spike rate asymmetry accurately lateralize the SOZ, '...
    'that relatively brief periods are needed to calculate spike rates, '...
    'and that the results are valid in an external dataset from a separate hospital.</p>']);


%% Add subtitles
annotation('textbox',[0 0.9 0.1 0.1],'String','A','LineStyle','none','fontsize',20)
annotation('textbox',[0.34 0.9 0.1 0.1],'String','B','LineStyle','none','fontsize',20)
annotation('textbox',[0.68 0.9 0.1 0.1],'String','C','LineStyle','none','fontsize',20)
annotation('textbox',[0 0.56 0.1 0.1],'String','D','LineStyle','none','fontsize',20)
annotation('textbox',[0.34 0.56 0.1 0.1],'String','E','LineStyle','none','fontsize',20)
annotation('textbox',[0.68 0.56 0.1 0.1],'String','F','LineStyle','none','fontsize',20)
annotation('textbox',[0 0.22 0.1 0.1],'String','G','LineStyle','none','fontsize',20)
annotation('textbox',[0.34 0.22 0.1 0.1],'String','H','LineStyle','none','fontsize',20)
annotation('textbox',[0.68 0.22 0.1 0.1],'String','I','LineStyle','none','fontsize',20)





print(gcf,[plot_folder,'Fig5'],'-dpng')

end