function model_plots(rm_wake,rm_non_temporal)

%% Parameters
pca_all_perc = 95;
pca_spikes_perc = 95;

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
set(gcf,'position',[1 1 1440 1000])
tiledlayout(2,3,"TileSpacing",'tight','padding','tight')

%% Run the lr_mt to extract features
if rm_wake == 1
    [T,features] =  lr_mt(3);
elseif rm_wake == 2
    [T,features] =  lr_mt(1);
else
    [T,features] =  lr_mt;
end
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];

%% ROC for L from R+BL for all features, just spikes, "binary spikes"
just_spikes = 0;% all patients (spikes for now for speed)
left = classifier_wrapper(T,features,pca_all_perc,1,just_spikes,rm_non_temporal,[]);
right = classifier_wrapper(T,features,pca_all_perc,2,just_spikes,rm_non_temporal,[]);

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
title({'Model performance (all features)'})
set(gca,'fontsize',20)

%% Feature importance
left_counts = left.counts;
left_unique = left.unique_features;
right_counts = right.counts;
right_unique = right.unique_features;

[~,Ileft] = sort(left_counts,'descend');
[~,Iright] = sort(right_counts,'descend');

% Plot
if 0
nexttile
bar(1:20,left_counts(Ileft(1:20)),'white')
text(1:20,left_counts(Ileft(1:20)),left_unique(Ileft(1:20)),...
    'horizontalalignment','right','Rotation',90,'fontsize',15,...
    'verticalalignment','middle')
xticklabels([])
xlabel('Features')
ylabel('Occurrences in top-ranked features')
title('Best features to identify left-sided SOZ')
set(gca,'fontsize',20)

end



%% ROC curve for multi-feature LR classifier to separate L from R+BL, also R from L+BL
% Do the classifier for L from R+BL
just_spikes = 1; % Just spikes
lefts = classifier_wrapper(T,features,pca_spikes_perc,1,just_spikes,rm_non_temporal,[]);
rights = classifier_wrapper(T,features,pca_spikes_perc,2,just_spikes,rm_non_temporal,[]);

% Get ROC stats
[XL,YL,~,AUCL] = perfcurve(lefts.class,lefts.scores,lefts.pos_class);
[XR,YR,~,AUCR] = perfcurve(rights.class,rights.scores,rights.pos_class);

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

%% Dumb spikes
just_spikes = 2; % Just spikes
leftd = classifier_wrapper(T,features,pca_spikes_perc,1,just_spikes,rm_non_temporal,[]);
rightd = classifier_wrapper(T,features,pca_spikes_perc,2,just_spikes,rm_non_temporal,[]);

% Get ROC stats
[XL,YL,~,AUCL] = perfcurve(leftd.class,leftd.scores,leftd.pos_class);
[XR,YR,~,AUCR] = perfcurve(rightd.class,rightd.scores,rightd.pos_class);

if 1
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
title({'Model performance (binary spikes)'})
set(gca,'fontsize',20)
end

%% Decide what method to use for further analyses
outl = lefts; % spikes model
outr = rights;

%% Compare model performance by hospital and by epilepsy localization
nb = 1e3;
%{
ADD BOOTSTRAPPING EVENTUALLY TO GET CONFIDENCE INTERVALS ON AUCS
%}
assert(isequal(outl.names,outr.names))



% Performance by hospital
hup = contains(outl.names,'HUP');
musc = contains(outl.names,'MP');
[~,~,~,AUCLH] = perfcurve(outl.class(hup),outl.scores(hup),outl.pos_class);
[~,~,~,AUCRH] = perfcurve(outr.class(hup),outr.scores(hup),outr.pos_class);
boot_lh = bootstrap_aucs(outl.class(hup),outl.scores(hup),outl.pos_class,nb);
boot_rh = bootstrap_aucs(outr.class(hup),outr.scores(hup),outr.pos_class,nb);
lh_ci = prctile(boot_lh,[2.5,97.5]);
rh_ci = prctile(boot_rh,[2.5,97.5]);


[~,~,~,AUCLM] = perfcurve(outl.class(musc),outl.scores(musc),outl.pos_class);
[~,~,~,AUCRM] = perfcurve(outr.class(musc),outr.scores(musc),outr.pos_class);
boot_lm = bootstrap_aucs(outl.class(musc),outl.scores(musc),outl.pos_class,nb);
boot_rm = bootstrap_aucs(outr.class(musc),outr.scores(musc),outr.pos_class,nb);
lm_ci = prctile(boot_lm,[2.5,97.5]);
rm_ci = prctile(boot_rm,[2.5,97.5]);

if ~rm_non_temporal
% Performance by SOZ localization
tle = contains(T.soz_locs,'temporal');
etle = contains(T.soz_locs,'other cortex') | contains(T.soz_locs,'frontal') | ...
    contains(T.soz_locs,'diffuse') | contains(T.soz_locs,'multifocal');

% Find the right patients
[~,~,~,AUCLT] = perfcurve(outl.class(tle),outl.scores(tle),outl.pos_class);
[~,~,~,AUCRT] = perfcurve(outr.class(tle),outr.scores(tle),outr.pos_class);
boot_lt = bootstrap_aucs(outl.class(tle),outl.scores(tle),outl.pos_class,nb);
boot_rt = bootstrap_aucs(outr.class(tle),outr.scores(tle),outr.pos_class,nb);
lt_ci = prctile(boot_lt,[2.5,97.5]);
rt_ci = prctile(boot_rt,[2.5,97.5]);

[~,~,~,AUCLE] = perfcurve(outl.class(etle),outl.scores(etle),outl.pos_class);
[~,~,~,AUCRE] = perfcurve(outr.class(etle),outr.scores(etle),outr.pos_class);
boot_le = bootstrap_aucs(outl.class(etle),outl.scores(etle),outl.pos_class,nb);
boot_re = bootstrap_aucs(outr.class(etle),outr.scores(etle),outr.pos_class,nb);
le_ci = prctile(boot_le,[2.5,97.5]);
re_ci = prctile(boot_re,[2.5,97.5]);
end

% plot
cols = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980];
nexttile
startg = 0.8;
smallg = 0.2;
mediumg = 0.4;
bigg = 0.8;
lp=errorbar(startg,AUCLH,AUCLH-lh_ci(1),lh_ci(2)-AUCLH,...
    'o','color',cols(1,:),'markersize',15,'linewidth',2);
hold on
rp = errorbar(startg+smallg,AUCRH,AUCRH-rh_ci(1),rh_ci(2)-AUCRH,...
    'o','color',cols(2,:),'markersize',15,'linewidth',2);

errorbar(startg+smallg+mediumg,AUCLM,AUCLM-lm_ci(1),lm_ci(2)-AUCLM,...
    'o','color',cols(1,:),'markersize',15,'linewidth',2)
hold on
errorbar(startg+2*smallg+mediumg,AUCRM,AUCRM-rm_ci(1),rm_ci(2)-AUCRM,...
    'o','color',cols(2,:),'markersize',15,'linewidth',2)

if ~rm_non_temporal
    errorbar(startg+2*smallg+mediumg+bigg,AUCLT,AUCLT-lt_ci(1),lt_ci(2)-AUCLT,...
        'o','color',cols(1,:),'markersize',15,'linewidth',2)
    hold on
    errorbar(startg+3*smallg+mediumg+bigg,AUCRT,AUCRT-rt_ci(1),rt_ci(2)-AUCRT,...
        'o','color',cols(2,:),'markersize',15,'linewidth',2)
    
    errorbar(startg+3*smallg+2*mediumg+bigg,AUCLE,AUCLE-le_ci(1),le_ci(2)-AUCLE,...
        'o','color',cols(1,:),'markersize',15,'linewidth',2)
    hold on
    errorbar(startg+4*smallg+2*mediumg+bigg,AUCRE,AUCRE-re_ci(1),re_ci(2)-AUCRE,...
        'o','color',cols(2,:),'markersize',15,'linewidth',2)
    xticks([(startg+startg+smallg)/2,(startg+smallg+mediumg+startg+2*smallg+mediumg)/2,...
        (startg+2*smallg+mediumg+bigg+startg+3*smallg+mediumg+bigg)/2,...
        (startg+3*smallg+2*mediumg+bigg+startg+4*smallg+2*mediumg+bigg)/2])
    xticklabels({'HUP','MUSC','TLE','ETLE'})
    xlim([startg-smallg startg+4*smallg+2*mediumg+bigg+smallg])
else
    xticks([(startg+startg+smallg)/2,(startg+smallg+mediumg+startg+2*smallg+mediumg)/2])
    xticklabels({'HUP','MUSC'})
    xlim([startg-smallg startg+2*smallg+mediumg+smallg])

end


ylabel('AUC')
leg_p = legend([lp,rp],{sprintf('Left vs right/bilateral'),...
    sprintf('Right vs left/bilateral')},'fontsize',20,...
    'location','south');
title('Model performance by clinical factors')
set(gca,'fontsize',20)


%% Validation
%{
train = contains(T.names,'HUP');
test = contains(T.names,'MP');
just_spikes = 1; % ALl features
rm_non_temporal = 0; % All patients
left = validation_classifier_wrapper(T,train,test,features,pca_perc,1,just_spikes,rm_non_temporal);
right = validation_classifier_wrapper(T,train,test,features,pca_perc,2,just_spikes,rm_non_temporal);
bilateral = validation_classifier_wrapper(T,train,test,features,pca_perc,3,just_spikes,rm_non_temporal);

% Get ROC stats
[XL,YL,~,AUCL] = perfcurve(left.class,left.scores,left.pos_class);
[XR,YR,~,AUCR] = perfcurve(right.class,right.scores,right.pos_class);so
[XB,YB,~,AUCB] = perfcurve(bilateral.class,bilateral.scores,bilateral.pos_class);

% Plot
nexttile
ll = plot(XL,YL,'linewidth',2);
hold on
lr = plot(XR,YR,':','linewidth',2);
b = plot(XB,YB,':','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr,b],{sprintf('Left vs right/bilateral: AUC = %1.2f',AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',AUCR),...
    sprintf('Bilateral vs left/right: AUC = %1.2f',AUCB)},'fontsize',20,...
    'location','southoutside')
title({'Model performance by laterality','Spike features, MUSC test dataset'})
set(gca,'fontsize',20)
%}


%

%% Three-way classifier
%{
out = classifier_wrapper(T,features,pca_perc,0,1,0);

% confusion matrix
nexttile
C = out.C;
classes = out.unique_classes;
nclasses = length(classes);
new_numbers = map_numbers_onto_range(C,[1 0]);
Ccolor = cat(3,ones(nclasses,nclasses,1),repmat(new_numbers,1,1,2));
D = diag(new_numbers);
Dcolor = [repmat(D,1,2),ones(length(D),1)];
Ccolor(logical(repmat(eye(nclasses,nclasses),1,1,3))) = Dcolor;
imagesc(Ccolor)
accuracy = sum(diag(C))/sum(C(:));
recall = nan(nclasses,1);
for i = 1:nclasses
    tp = C(i,i);
    fn = sum(C(i,~ismember(1:nclasses,i))); 
    recall(i) = tp/(tp+fn); % tp is number correctly predicted to be in class, tp + fn is everyone in the class
end
balanced_accuracy = mean(recall);

xticks(1:nclasses)
xticklabels((classes))
yticks(1:nclasses)
yticklabels((classes))
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
%}
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
title(sprintf('Accuracy:%1.1f%%\nBalanced accuracy: %1.1f%%',...
    accuracy*100,balanced_accuracy*100))
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



print(gcf,[plot_folder,'Fig5'],'-dpng')

end