function trainedClassifier = pca_then_fs_classifier(trainingData,method,features,respVar,pc_perc,ncycles,numFeaturesToKeep)

%% Data prep
inputTable = trainingData;
predictorNames = features;
predictors = inputTable(:, predictorNames);
response = inputTable.(respVar);
classNames = unique(response);

%% PCA - figure out how to normalize
% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
numericPredictors = predictors;
numericPredictors = table2array(varfun(@double, numericPredictors));
% 'inf' values have to be treated as missing data for PCA.
numericPredictors(isinf(numericPredictors)) = NaN;

% Do PCA
w = 1./var(numericPredictors,[],1,"omitnan");
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors,'centered',true,'VariableWeights',w);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = pc_perc/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep))];


%% 5-fold cross validation to select features
KFolds = 5;
cvp = cvpartition(response, 'KFold', KFolds);
all_included_predictors = cell(KFolds,1);

for fold = 1:KFolds

    temp_predictors = predictors(cvp.training(fold), :);
    temp_response = response(cvp.training(fold), :);

    % Feature Ranking and Selection
    % Replace Inf/-Inf values with NaN to prepare data for normalization
    temp_predictors = standardizeMissing(temp_predictors, {Inf, -Inf});
    predictorMatrix = array2table((table2array(temp_predictors)-nanmean(table2array(temp_predictors),1))./nanstd(table2array(temp_predictors),[],1));
    newPredictorMatrix = zeros(size(predictorMatrix));
    for i = 1:size(predictorMatrix, 2)
        newPredictorMatrix(:,i) = predictorMatrix{:,i}; 
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
    tempIncludedPredictorNames = temp_predictors.Properties.VariableNames(featureIndex(1:numFeaturesToKeep));
    all_included_predictors{fold} = tempIncludedPredictorNames;
end
all_included_predictors = horzcat(all_included_predictors{:});

% Take the N with the most votes
a=unique(all_included_predictors,'stable');
b=cellfun(@(x) sum(ismember(all_included_predictors,x)),a,'un',0); % get counts for each
counts = cell2mat(b);

% Make sure counts are sorted
[~,I] = sort(counts,'descend');
top_predictors = a(I(1:numFeaturesToKeep));
includedPredictorNames = top_predictors;

% Assign these as the new predictors
predictors = predictors(:,includedPredictorNames);




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
        %{
        template = templateTree(...
            'MaxNumSplits', 100, ...
            'NumVariablesToSample', 'all');
        %}
        template = templateTree('PredictorSelection','interaction-curvature');
        
        classifier = fitcensemble(...
            predictors, ...
            response, ...
            'Method', 'bag', ...
            'Learners', template, ...
            'ClassNames', classNames,...
            'NumLearningCycles',ncycles);
    case 'boost'
        
        template = templateTree(...
            'MaxNumSplits', 10, ...
            'NumVariablesToSample', 'all',...
            'PredictorSelection','interaction-curvature');
        
        classifier = fitcensemble(...
            predictors, ...
            response, ...
            'Method', 'AdaBoostM2', ...
            'Learners', template, ...
            'ClassNames', classNames,...
            'NumLearningCycles',ncycles);

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
    case 'lr' % only for binary
        classifier = fitclinear(predictors,response);

    case 'naiveBayes'
        distributionNames =  repmat({'Normal'}, 1, size(predictors,2));
        classifier = fitcnb(...
            predictors, ...
            response, ...
            'DistributionNames', distributionNames, ...
            'ClassNames', classNames);
       
end

%% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
featureSelectionFcn = @(x) x(:,includedPredictorNames);
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x)) - pcaCenters) .* w * pcaCoefficients)];
oldPredictFcn = @(x) predict(classifier, x);
predictFcn = @(x) oldPredictFcn(featureSelectionFcn(pcaTransformationFcn(predictorExtractionFcn(x))));

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



