function trainedClassifier = class_feature_select(trainingData,method,features,respVar,pc_perc,ncycles)

%% Data prep
inputTable = trainingData;
predictorNames = features;
predictors = inputTable(:, predictorNames);
response = inputTable.(respVar);
isCategoricalPredictor = repmat(false,1,length(predictorNames));
classNames = unique(response);

%nca = fscnca(table2array(predictors),response,'FitMethod','exact');

%% 5-fold cross validation to select optimal features
KFolds = 5;
cvp = cvpartition(response, 'KFold', KFolds);
numvalidsets = cvp.NumTestSets;
nFeatureVals = 5:30;
lossVals = zeros(length(nFeatureVals),numvalidsets);
all_included_predictors = cell(length(nFeatureVals),numvalidsets);

for k = 1:length(nFeatureVals)

    for fold = 1:KFolds
    
        temp_predictors = predictors(cvp.training(fold), :);
        temp_response = response(cvp.training(fold), :);
    
        % Feature Ranking and Selection
        % Replace Inf/-Inf values with NaN to prepare data for normalization
        temp_predictors = standardizeMissing(temp_predictors, {Inf, -Inf});
        % Normalize data for feature ranking
        predictorMatrix = normalize(temp_predictors, "DataVariable", ~isCategoricalPredictor);
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
        pvalues = nan(size(predictorMatrix, 2),1);
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
        all_included_predictors{k,fold} = tempIncludedPredictorNames;

        % Do a decision tree on training data
        classifier = fitctree(...
            temp_selected_predictors, ...
            temp_response, ...
            'SplitCriterion', 'gdi', ...
            'MaxNumSplits', 100, ...
            'Surrogate', 'off', ...
            'ClassNames', classNames);

        % test on testing data
        test_response = response(cvp.testing(fold), :);


    end
    
    

    

end


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
    numericPredictors);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = pc_perc/100;
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



