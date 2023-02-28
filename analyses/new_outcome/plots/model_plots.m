function model_plots

%% Parameters
pca_perc = 90; %90 works

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


%% Initialize figure
figure
set(gcf,'position',[1 1 1400 450])
tiledlayout(1,3,"TileSpacing",'tight','padding','tight')

%% Run the lr_mt to extract features
[T,features] =  lr_mt;

%% A: ROC curve for multi-feature LR classifier to separate L from R+BL, also R from L+BL
% Do the classifier for L from R+BL
just_spikes = 0; % ALl features
rm_non_temporal = 0; % All patients
left = classifier_wrapper(T,features,pca_perc,1,just_spikes,rm_non_temporal);
right = classifier_wrapper(T,features,pca_perc,2,just_spikes,rm_non_temporal);

% Get ROC stats
[XL,YL,~,AUCL] = perfcurve(left.class,left.scores,left.pos_class);
[XR,YR,~,AUCR] = perfcurve(right.class,right.scores,right.pos_class);

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
title('Model performance by laterality')
set(gca,'fontsize',20)

%% B: ROC for L from R+BL for all features, just spikes, "binary spikes"
combine_br = 1; % L vs R+BL
rm_non_temporal = 0; % all patients
all = classifier_wrapper(T,features,pca_perc,combine_br,0,rm_non_temporal);
spikes = classifier_wrapper(T,features,pca_perc,combine_br,1,rm_non_temporal);
bin_spikes = classifier_wrapper(T,features,pca_perc,combine_br,2,rm_non_temporal);

% Get ROC stats
[XA,YA,~,AUCA] = perfcurve(all.class,all.scores,all.pos_class);
[XS,YS,~,AUCS] = perfcurve(spikes.class,spikes.scores,spikes.pos_class);
[XB,YB,~,AUCB] = perfcurve(bin_spikes.class,bin_spikes.scores,bin_spikes.pos_class);

% Plot
nexttile
la = plot(XA,YA,'linewidth',2);
hold on
ls = plot(XS,YS,':','linewidth',2);
lb = plot(XB,YB,'-.','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([la,ls,lb],{sprintf('All features: AUC = %1.2f',AUCA),...
    sprintf('Spike features: AUC = %1.2f',AUCS),sprintf('Spike features binarized: AUC = %1.2f',AUCB)},'fontsize',20,...
    'location','southeast')
title('Model performance by features')
set(gca,'fontsize',20)


%% C: ROC for L from R+BL for temporal vs ETLE - no don't do this
%{
% There aren't enough patients who are ETLE for this to be a meaningful
% analysis. AUC is like 0.25 but for the number of patients this is not
% significant.
just_spikes = 0;
combine_br = 1;
temporal = classifier_wrapper(T,features,pca_perc,combine_br,just_spikes,1);
extra = classifier_wrapper(T,features,pca_perc,combine_br,just_spikes,2);

% Get ROC stats
[XT,YT,~,AUCT] = perfcurve(temporal.class,temporal.scores,temporal.pos_class);
[XE,YE,~,AUCE] = perfcurve(extra.class,extra.scores,extra.pos_class);

% investigating ETLE
n1 = sum(strcmp(extra.class,'br'));
n2 = sum(strcmp(extra.class,'left'));
U = AUCE*n1*n2;
mu = n1*n2/2;
sigmau = sqrt(n1*n2*(n1+n2+1)/12);
z = (U-mu)/sigmau;
p = normcdf(z);

% Plot
nexttile
lt = plot(XT,YT,':','linewidth',2);
hold on
le = plot(XE,YE,'linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([lt,le],{sprintf('Temporal (N = %d): AUC = %1.2f',temporal.npts,AUCT),...
    sprintf('Extra-temporal (N = %d): AUC = %1.2f',extra.npts,AUCE)},'fontsize',20,...
    'location','southeast')
title('Model performance by SOZ localization')
set(gca,'fontsize',20)
%}

%% D: Confusion matrix for threshold 0.5
C = all.C;
classes = all.unique_classes;
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
title(sprintf('Accuracy: %1.1f%%, balanced accuracy: %1.1f%%',...
    accuracy*100,balanced_accuracy*100))
set(gca,'fontsize',20)

print(gcf,[plot_folder,'Fig5'],'-dpng')

end