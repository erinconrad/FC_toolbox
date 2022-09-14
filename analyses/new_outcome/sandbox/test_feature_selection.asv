function test_feature_selection


%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Basic parameters
N = 1e2;
nobservations = 50;
perc_train = 2/3;
npredictors = 0;
model = @(x,y,z) fitglm(x,'ResponseVar',y,'PredictorVars',z,'Distribution','binomial');

%% Make coin flip data (CV accuracy should be 50%)
T=makeCoinFlipData(nobservations,npredictors,'binary','binary');
response = 'response';
predictors = T.Properties.VariableNames(~strcmp(T.Properties.VariableNames,response));

%% First, show that CV error increases as you add parameters
auc_n_predictors = nan(npredictors,1);
for i = 1:npredictors
    [~,~,AUC] = more_general_train_test(T,N,perc_train,model,predictors(1:i),response);
    auc_n_predictors(i) = nanmean(AUC);
end


%% Now, show that you reduce cross validation error if you create BIAS by doing feature selection FIRST
% Do feature selection
idx = fscchi2(T,response);

auc_best_n_predictors = nan(npredictors,1);
for i = 1:npredictors
    best_predictors = predictors(idx(1:i));
    [~,~,AUC] = more_general_train_test(T,N,perc_train,model,best_predictors,response);
    auc_best_n_predictors(i) = nanmean(AUC);
end

%% Next, show that if you instead do feature selection within each fold, you avoid bias
auc_best_n_predictors_within_fold = nan(npredictors,1);
for i = 1:npredictors
    [~,~,AUC] = train_test_chi2fs(T,N,perc_train,model,predictors,response,i);
    auc_best_n_predictors_within_fold(i) = nanmean(AUC);
end



figure
tiledlayout(1,3)
nexttile
plot(1:npredictors,auc_n_predictors,'o','linewidth',2,'markersize',10)
hold on
plot(xlim,[mean(auc_n_predictors) mean(auc_n_predictors)],'--','linewidth',2)
xlabel('Number of predictors');ylabel('Cross validation accuracy')
ylim([0 1])
title('All predictors')
set(gca,'fontsize',15)

nexttile
plot(1:npredictors,auc_best_n_predictors,'o','linewidth',2,'markersize',10)
hold on
plot(xlim,[mean(auc_best_n_predictors) mean(auc_best_n_predictors)],'--','linewidth',2)
xlabel('Number of predictors');ylabel('Cross validation accuracy')
ylim([0 1])
title('Best predictors')
set(gca,'fontsize',15)

nexttile
plot(1:npredictors,auc_best_n_predictors_within_fold,'o','linewidth',2,'markersize',10)
hold on
plot(xlim,[mean(auc_best_n_predictors_within_fold) mean(auc_best_n_predictors_within_fold)],'--','linewidth',2)
xlabel('Number of predictors');ylabel('Cross validation accuracy')
ylim([0 1])
title('Best predictors (selection in fold)')
set(gca,'fontsize',15)
%}

end