function test_external_validation

% This code performs both leave-one-out cross validation on the HUP
% training data AND trains the model on all the HUP data and tests it on
% the MUSC data. It does this for a model using all features as well as a
% model using only spike rates.

%% Parameters
for_will = 0;
pca_perc = 95; % the percent variance to explain for pca
rm_non_temporal = 1; % remove patients who are not temporal
rm_wake = 1; % don't include wake segments
rm_bad_spikes = 0;
which_ref = 'car';

%% Get file locs

locations = fc_toolbox_locs;
data_folder = locations.el_data_folder;
fmri_folder = [data_folder,'fmri_data/'];
plot_folder = locations.el_plots_folder;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load the file containing intermediate data
inter_folder = data_folder;
mt_data = load([inter_folder,'mt_out_epilepsy_laterality.mat']);
mt_data = mt_data.out;

%% RID table
rid_table = readtable([inter_folder,'all_rids.csv']);

%% fmri locs
file_path = fmri_folder;
csv_path = [file_path,'out_csvs/'];
fT = readtable([file_path,'df.csv']);

%% Run the lr_mt to extract AI features
if rm_wake == 1
    [T,features] =  lr_mt(mt_data,3,rm_bad_spikes); % the 3 refers to only looking at sleep

else
    error('why are you doing this?')
end



% Remove those without a response (soz_lats is the response variable)
empty_class = cellfun(@isempty,T.soz_lats);
T(empty_class,:) = [];

% As a test, remove bilateral
if 0
    bilateral = strcmp(T.soz_lats,'bilateral');
    T(bilateral,:) = [];
end

%% Get fmri AI and add it to table
T = add_fmri_info_table(T,fT,rid_table,csv_path);

%% Establish HUP and MUSC as training and testing, respectively
train  = contains(T.names,'HUP');
test  = contains(T.names,'MP');

%% try multiclass model on spikes
% oof doesn't work
%{
Ttrain = T(train,:);
just_spikes = 1;
sint = classifier_wrapper(Ttrain,features,pca_perc,0,just_spikes,rm_non_temporal,[]); % 0 means multiclass
sext = validation_classifier_wrapper(T,train,test,features,pca_perc,0,just_spikes,rm_non_temporal); % 1 means left
%}

%% Do the LOO cross validation on the HUP data - FULL model
Ttrain = T(train,:);
just_spikes = 0; % not just spikes, full model
lefta_int = classifier_wrapper(Ttrain,features,pca_perc,1,just_spikes,rm_non_temporal,[],which_ref); % 1 means left
righta_int = classifier_wrapper(Ttrain,features,pca_perc,2,just_spikes,rm_non_temporal,[],which_ref); % 2 means right

% Get ROC stats
[lefta_int.XL,lefta_int.YL,~,lefta_int.AUCL] = perfcurve(lefta_int.class,lefta_int.scores,lefta_int.pos_class);
[righta_int.XR,righta_int.YR,~,righta_int.AUCR] = perfcurve(righta_int.class,righta_int.scores,righta_int.pos_class);

%% Train on the HUP data, test on MUSC - FULL model
just_spikes = 0; % not just spikes
lefta_ext = validation_classifier_wrapper(T,train,test,features,pca_perc,1,just_spikes,rm_non_temporal,which_ref); % 1 means left
righta_ext = validation_classifier_wrapper(T,train,test,features,pca_perc,2,just_spikes,rm_non_temporal,which_ref); % 2 means right

% Get ROC stats
[lefta_ext.XL,lefta_ext.YL,~,lefta_ext.AUCL] = perfcurve(lefta_ext.class,lefta_ext.scores,lefta_ext.pos_class);
[righta_ext.XR,righta_ext.YR,~,righta_ext.AUCR] = perfcurve(righta_ext.class,righta_ext.scores,righta_ext.pos_class);


%% Do the LOO cross validation on the HUP data - spike only model
Ttrain = T(train,:);
just_spikes = 1; % only spike feature
lefts_int = classifier_wrapper(Ttrain,features,pca_perc,1,just_spikes,rm_non_temporal,[],which_ref); % 1 means left
rights_int = classifier_wrapper(Ttrain,features,pca_perc,2,just_spikes,rm_non_temporal,[],which_ref); % 2 means right

% Get ROC stats
[lefts_int.XL,lefts_int.YL,~,lefts_int.AUCL] = perfcurve(lefts_int.class,lefts_int.scores,lefts_int.pos_class);
[rights_int.XR,rights_int.YR,~,rights_int.AUCR] = perfcurve(rights_int.class,rights_int.scores,rights_int.pos_class);


%% Train on the HUP data, test on MUSC - spike only model
fprintf('\nDoing main models...');
just_spikes = 1; 
lefts_ext = validation_classifier_wrapper(T,train,test,features,pca_perc,1,just_spikes,rm_non_temporal,which_ref); % 1 means left
rights_ext = validation_classifier_wrapper(T,train,test,features,pca_perc,2,just_spikes,rm_non_temporal,which_ref); % 2 means right

% Get ROC stats
[lefts_ext.XL,lefts_ext.YL,~,lefts_ext.AUCL] = perfcurve(lefts_ext.class,lefts_ext.scores,lefts_ext.pos_class);
[rights_ext.XR,rights_ext.YR,~,rights_ext.AUCR] = perfcurve(rights_ext.class,rights_ext.scores,rights_ext.pos_class);


%% Show the results
figure
tiledlayout(2,2)

% Full model - internal
nexttile
ll = plot(lefta_int.XL,lefta_int.YL,'linewidth',2);
hold on
lr = plot(righta_int.XR,righta_int.YR,':','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',lefta_int.AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',righta_int.AUCR)},'fontsize',20,...
    'location','southeast')
title({'Full model - LOO CV HUP'})
set(gca,'fontsize',20)

% Full model - external
nexttile
ll = plot(lefta_ext.XL,lefta_ext.YL,'linewidth',2);
hold on
lr = plot(righta_ext.XR,righta_ext.YR,':','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',lefta_ext.AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',righta_ext.AUCR)},'fontsize',20,...
    'location','southeast')
title({'Full model - External validation'})
set(gca,'fontsize',20)

% Spike model - internal
nexttile
ll = plot(lefts_int.XL,lefts_int.YL,'linewidth',2);
hold on
lr = plot(rights_int.XR,rights_int.YR,':','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',lefts_int.AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',rights_int.AUCR)},'fontsize',20,...
    'location','southeast')
title({'Spike model - LOO CV HUP'})
set(gca,'fontsize',20)

% Spike model - external
nexttile
ll = plot(lefts_ext.XL,lefts_ext.YL,'linewidth',2);
hold on
lr = plot(rights_ext.XR,rights_ext.YR,':','linewidth',2);
plot([0 1],[0 1],'k--','linewidth',2)
xlabel('False positive rate')
ylabel('True positive rate')
legend([ll,lr],{sprintf('Left vs right/bilateral: AUC = %1.2f',lefts_ext.AUCL),...
    sprintf('Right vs left/bilateral: AUC = %1.2f',rights_ext.AUCR)},'fontsize',20,...
    'location','southeast')
title({'Spike model - External validation'})
set(gca,'fontsize',20)

%% double check some things - these checks work
% basically, I confirmed that I can use a simple logistic function that
% takes the coefficients from the model and re-derive the scores, and then
% if I apply a cutoff of 0.5 it derives the predictions.
% get info for left model
features_left_test = T{test,sprintf('spikes %s sleep',which_ref)};
classNames = lefts_ext.unique_classes;
coefs = lefts_ext.tc.coef;

% Try to solve it again using a simple calculator
[preds,scores] = simple_prediction_calculator(features_left_test,classNames,coefs);

% confirm that I can re-derive the predicted scores and classes
assert(isequal(preds,lefts_ext.all_pred))
assert(sum(abs(scores-lefts_ext.scores)>1e-3)==0)

% same for right
features_right_test = T{test,sprintf('spikes %s sleep',which_ref)};
classNames = rights_ext.unique_classes;
coefs = rights_ext.tc.coef;

% Try to solve it again using a simple calculator
[preds,scores] = simple_prediction_calculator(features_right_test,classNames,coefs);

% confirm that I can re-derive the predicted scores and classes
assert(isequal(preds,rights_ext.all_pred))
assert(sum(abs(scores-rights_ext.scores)>1e-3)==0)
end