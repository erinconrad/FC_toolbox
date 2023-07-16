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
out = load([plot_folder,'models.mat']);
out = out.out;

%% Initialize results file
fname = [plot_folder,'results.html'];
fid = fopen(fname,'a');
fprintf(fid,'<p><br><b>A classifier incorporating interictal EEG features predicts epilepsy laterality</b></br>');



%% Initialize figure
figure
set(gcf,'position',[1 1 1440 1000])
tiledlayout(2,3,"TileSpacing",'tight','padding','tight')

%% ROC for all features
% Get ROC curve stuff
left = out.full_model.left;
right = out.full_model.right;
XL = left.XL; YL = left.YL; AUCL = left.AUCL;
XR = right.XR; YR = right.YR; AUCR = right.AUCR;

% Plot
nexttile
ll = plot(XL,YL,'linewidth',2);
hold on
lr = plot(XR,YR,':','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',AUCR)},'fontsize',20,...
    'location','southeast')
title({'Model performance (all features)'})
set(gca,'fontsize',20)

fprintf(fid,['We trained two models incorporating all interictal features to predict '...
    'SOZ laterality, with one model to predict left versus right/bilateral '...
    'onset and one model to predict right versus left/bilateral onset. The '...
    'area under the curve (AUC) of the receiver operator characteristic (ROC) curve '...
    'was %1.2f for the left-sided classifier and %1.2f for the right-sided classifier '...
    '(Fig. 5A).'],AUCL,AUCR);


%% ROC curve just spikes
% Get ROC stuff
lefts = out.spike_model.left;
rights = out.spike_model.right;
XL = lefts.XL; YL = lefts.YL; AUCL = lefts.AUCL;
XR = rights.XR; YR = rights.YR; AUCR = rights.AUCR;

% Plot
nexttile
ll = plot(XL,YL,'linewidth',2);
hold on
lr = plot(XR,YR,':','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',AUCR)},'fontsize',20,...
    'location','southeast')
title({'Model performance (spikes only)'})
set(gca,'fontsize',20)

fprintf(fid,['Next, we examined how a model trained using just average spikes rates '...
    '(bipolar reference) would perform. We chose this single feature because this is '...
    'a feature that could be manually derived using current clinical approaches (e.g., counting '...
    'spikes on the left and on the right). The '...
    'area under the curve (AUC) of the receiver operator characteristic (ROC) curve '...
    'was %1.2f for the left-sided classifier and %1.2f for the right-sided classifier '...
    '(Fig. 5A).'],AUCL,AUCR);

%% ROC curve binary spikes
% Get ROC stuff
leftd = out.binary_spike_model.left;
rightd = out.binary_spike_model.right;
XL = leftd.XL; YL = leftd.YL; AUCL = leftd.AUCL;
XR = rightd.XR; YR = rightd.YR; AUCR = rightd.AUCR;


nexttile
ll = plot(XL,YL,'linewidth',2);
hold on
lr = plot(XR,YR,':','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',AUCR)},'fontsize',20,...
    'location','southeast')
title({'Model performance (binary spikes)'})
set(gca,'fontsize',20)

%% Decide what method to use for further analyses
outl = lefts; % spikes model
outr = rights;


%% Confusion matrix for threshold 0.5 for left (spikes only)
C = outl.C;
classes = outl.unique_classes;
nclasses = length(classes);

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
pretty_name = strrep(pretty_name,'br','Right/bilateral');
xticks(1:nclasses)
xticklabels((pretty_name))
yticks(1:nclasses)
yticklabels((pretty_name))
xlabel('Predicted')
ylabel('True')
hold on
for i = 1:nclasses
    for j = 1:nclasses
        text(i,j,sprintf('%d',C(j,i)),'horizontalalignment','center','fontsize',20)
    end
end
title(sprintf('Example lateralization using threshold AI = 0.4\nAccuracy: %1.1f (post-hoc, not cross-validated)',accuracy*100))
set(gca,'fontsize',20)


%% Confusion matrix for threshold 0.5 for right (spikes only)
C = outr.C;
classes = outr.unique_classes;
nclasses = length(classes);

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
pretty_name = strrep(pretty_name,'right','Right');
pretty_name = strrep(pretty_name,'bl','Left/bilateral');
xticks(1:nclasses)
xticklabels((pretty_name))
yticks(1:nclasses)
yticklabels((pretty_name))
xlabel('Predicted')
ylabel('True')
hold on
for i = 1:nclasses
    for j = 1:nclasses
        text(i,j,sprintf('%d',C(j,i)),'horizontalalignment','center','fontsize',20)
    end
end
title(sprintf('Accuracy: %1.1f%%\nBalanced accuracy: %1.1f%%',...
    accuracy*100,balanced_accuracy*100))
set(gca,'fontsize',20)

%% Subsampling plots
nexttile

sub = out.subsampling.data;
durations = out.subsampling.durations;
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

ylim([0.6 0.9])

legend([el,er],{'Left vs right/bilateral','Right vs left/bilateral'},'location','southeast')
xticks(1:ndurs)
xticklabels(arrayfun(@(x) sprintf('%d min',x),durations,'uniformoutput',false))
ylabel('Median (IQR) AUC')
title(sprintf('Model accuracy by duration'))
set(gca,'fontsize',20)


print(gcf,[plot_folder,'Fig5'],'-dpng')

end