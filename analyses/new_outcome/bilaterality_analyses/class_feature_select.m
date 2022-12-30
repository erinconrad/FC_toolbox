function trainedClassifier = class_feature_select(trainingData,method,features,respVar,ncycles)

%% Data prep
inputTable = trainingData;
predictorNames = features;
predictors = inputTable(:, predictorNames);
response = inputTable.(respVar);
isCategoricalPredictor = repmat(false,1,length(predictorNames));
classNames = unique(response);

%nca = fscnca(table2array(predictors),response,'FitMethod','exact');

%% 5-fold cross validation to select optimal 1) # of features, 2) features, 3) PCA 
KFolds = 5;
cvp = cvpartition(response, 'KFold', KFolds);
numvalidsets = cvp.NumTestSets;
nFeatureVals = 2:20;
pca_percs = 90;
lossVals = zeros(length(nFeatureVals),length(pca_percs),numvalidsets);
all_included_predictors = cell(length(nFeatureVals),length(pca_percs),numvalidsets);

% Loop over number of features
for k = 1:length(nFeatureVals)

    % loop over pca percent
    for j = 1:length(pca_percs)

        % loop over folds
        for fold = 1:KFolds
        
            % Generate fold
            temp_predictors = predictors(cvp.training(fold), :);
            temp_response = response(cvp.training(fold), :);
        
            %% Feature Ranking and Selection
            % Replace Inf/-Inf values with NaN to prepare data for normalization
            temp_predictors = standardizeMissing(temp_predictors, {Inf, -Inf});

            % Normalize data for feature ranking
            %predictorMatrix = normalize(temp_predictors, "DataVariable", ~isCategoricalPredictor);
             predictorMatrix = array2table((table2array(temp_predictors)-nanmean(table2array(temp_predictors),1))./nanstd(table2array(temp_predictors),[],1));

            newPredictorMatrix = zeros(size(predictorMatrix));
            for i = 1:size(predictorMatrix, 2)
                if isCategoricalPredictor(i)
                    newPredictorMatrix(:,i) = grp2idx(predictorMatrix{:,i});
                else
                    newPredictorMatrix(:,i) = predictorMatrix{:,i};
                end
            end
            predictorMatrix = newPredictorMatrix;
            responseVector = grp2idx(temp_response);
            
            % Rank features using Kruskal Wallis algorithm
            pValues = nan(size(predictorMatrix, 2),1);
            for i = 1:size(predictorMatrix, 2)
                pValues(i) = kruskalwallis(...
                    predictorMatrix(:,i), ...
                    responseVector, ...
                    'off');
            end
            [~,featureIndex] = sort(-log(pValues), 'descend');
            tempIncludedPredictorNames = temp_predictors.Properties.VariableNames(featureIndex(1:nFeatureVals(k)));
            temp_selected_predictors = temp_predictors(:,tempIncludedPredictorNames);
    
            % store the predictors
            all_included_predictors{k,j,fold} = tempIncludedPredictorNames;

            %% PCA            
            % Normalize the predictors prior to PCA
            normalizationFcn = @(x) (x-nanmean(x,1))./nanstd(x,[],1);
            normalizedPredictors = normalizationFcn(table2array(temp_selected_predictors));
            
            % Do PCA
            [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
                normalizedPredictors);
            % Keep enough components to explain the desired amount of variance.
            explainedVarianceToKeepAsFraction = pca_percs(j)/100;
            numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
            pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
            temp_post_pca_predictors = array2table(pcaScores(:,1:numComponentsToKeep));

    
            %% Do a decision tree on training data
            classifier = fitctree(...
                temp_post_pca_predictors, ...
                temp_response, ...
                'SplitCriterion', 'gdi', ...
                'MaxNumSplits', 100, ...
                'Surrogate', 'off', ...
                'ClassNames', classNames);
    
            %% test on testing data
            test_response = response(cvp.test(fold));
            test_predictors = predictors(cvp.test(fold),:);
            test_predictors_post_feature_selection = test_predictors(:,tempIncludedPredictorNames);
            test_predictors_post_pca = array2table((table2array(test_predictors_post_feature_selection) - pcaCenters) * pcaCoefficients);
            test_prediction = predict(classifier,test_predictors_post_pca);
            
            % define loss to be 1-accuracy
            temp_loss = sum(cellfun(@(x,y) ~strcmp(x,y),test_response,test_prediction))/length(test_response);
            lossVals(k,j,fold) = temp_loss;
    
        end
    
    
    end
    
end

%% Average the loss across folds
avg_loss = mean(lossVals,3);

if 0
    figure
    turn_nans_gray(avg_loss)
    colorbar
    xticks(1:length(pca_percs))
    yticks(1:length(nFeatureVals))
    xticklabels(pca_percs)
    yticklabels(nFeatureVals)
end

% Find the pca perc and num features corresponding to the lowest loss
[M,I] = min(avg_loss,[],'all');
[r,c] = ind2sub(size(avg_loss),I);
best_pca_perc = pca_percs(c);
best_nfeatures = nFeatureVals(r);

%% Now find the best features for this number of features
feature_set = horzcat(all_included_predictors{r,c,:});
% Take the N with the most votes
a=unique(feature_set,'stable');
b=cellfun(@(x) sum(ismember(feature_set,x)),a,'un',0); % get counts for each
counts = cell2mat(b);

% Make sure counts are sorted
[~,I] = sort(counts,'descend');
top_predictors = a(I(1:best_nfeatures));
includedPredictorNames = top_predictors;

% Assign these as the new predictors
predictors = predictors(:,includedPredictorNames);
isCategoricalPredictor = repmat(false,1,length(includedPredictorNames));

%% PCA - figure out how to normalize


% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 'inf' values have to be treated as missing data for PCA.
numericPredictors(isinf(numericPredictors)) = NaN;

% Normalize the predictors prior to PCA
normalizationFcn = @(x) (x-nanmean(x,1))./nanstd(x,[],1);
normalizedPredictors = normalizationFcn(numericPredictors);

if 0
    C = corr(numericPredictors);
    turn_nans_gray(C)
    colorbar
    xticklabels(includedPredictorNames)
    yticklabels(includedPredictorNames)
    set(gca,'fontsize',15)
end

% Do PCA
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    normalizedPredictors);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = best_pca_perc/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep)), predictors(:, isCategoricalPredictor)];

% Do some checks about what comprises the main components
if 0
    figure
    biplot(pcaCoefficients(:,1:2),'Scores',pcaScores(:,1:2),'VarLabels',includedPredictorNames)
end

%% Train a classifier
% This code specifies all the classifier options and trains the classifier.
switch method
    case 'tree'
        classifier = fitctree(...
            predictors, ...
            response, ...
            'SplitCriterion', 'gdi', ...
            'MaxNumSplits', 100, ...
            'Surrogate', 'off', ...
            'ClassNames', classNames);
    case 'knn'
        classifier = fitcknn(...
            predictors, ...
            response, ...
            'Distance', 'Euclidean', ...
            'Exponent', [], ...
            'NumNeighbors', 10, ...
            'DistanceWeight', 'Equal', ...
            'Standardize', true, ...
            'ClassNames', classNames);
    case 'bag'
        template = templateTree(...
            'MaxNumSplits', 100, ...
            'NumVariablesToSample', 'all');
        
        classifier = fitcensemble(...
            predictors, ...
            response, ...
            'Method', 'bag', ...
            'Learners', template, ...
            'ClassNames', classNames,...
            'NumLearningCycles',ncycles);
    case 'boost'
        template = templateTree(...
            'MaxNumSplits', 100, ...
            'NumVariablesToSample', 'all');
        
        classifier = fitcensemble(...
            predictors, ...
            response, ...
            'Method', 'AdaBoostM2', ...
            'Learners', template, ...
            'ClassNames', classNames,...
            'NumLearningCycles',100);
    case 'fancy_bag'
        template = templateTree(...
            'MaxNumSplits', 100, ...
            'NumVariablesToSample', 'all');
        
        classifier = fitcensemble(predictors,...
            response, ...
            'OptimizeHyperparameters','auto','Learners',template, ...
            'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus'));
    case 'svm'

        template = templateSVM(...
            'KernelFunction', 'linear', ...
            'PolynomialOrder', [], ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true);
        classifier = fitcecoc(...
            predictors, ...
            response, ...
            'Learners', template, ...
            'Coding', 'onevsone', ...
            'ClassNames', classNames);

    case 'naiveBayes'
        distributionNames =  repmat({'Normal'}, 1, numComponentsToKeep);
        classifier = fitcnb(...
            predictors, ...
            response, ...
            'DistributionNames', distributionNames, ...
            'ClassNames', classNames);
       
end

%% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
featureSelectionFcn = @(x) x(:,includedPredictorNames);
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
oldPredictFcn = @(x) predict(classifier, x);
predictFcn = @(x) oldPredictFcn(pcaTransformationFcn(featureSelectionFcn(predictorExtractionFcn(x))));

% Add additional fields to the result struct
trainedClassifier.predictFcn = predictFcn;
trainedClassifier.RequiredVariables = predictorNames;
trainedClassifier.PCACenters = pcaCenters;
trainedClassifier.PCACoefficients = pcaCoefficients;
trainedClassifier.classifier = classifier;
trainedClassifier.pcaTransformationFcn = pcaTransformationFcn;
trainedClassifier.featureSelectionFcn = featureSelectionFcn;
trainedClassifier.predictorExtractionFcn = predictorExtractionFcn;
trainedClassifier.oldPredictFcn = oldPredictFcn;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2022a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');



