function trainedClassifier = nca_classifier(trainingData,method,features,respVar,pca_perc,ncycles)

%% Data prep
inputTable = trainingData;
predictorNames = features;
predictors = inputTable(:, predictorNames);
response = inputTable.(respVar);
isCategoricalPredictor = repmat(false,1,length(predictorNames));
classNames = unique(response);


%% NCA for feature selection
% Get data into format needed for NCA
ytrain = response;
Xtrain = table2array(predictors);

% 5 fold cross validation
KFolds = 5;
cvp = cvpartition(ytrain, 'KFold', KFolds);
numvalidsets = cvp.NumTestSets;
n = length(ytrain);
lambdavals = linspace(0,4,10)/n;
lossvals = zeros(length(lambdavals),numvalidsets);

% Loop over lambdas
for i = 1:length(lambdavals)

    % Loop over folds
    for k = 1:numvalidsets

        % Get fold-level training and testing data
        X = Xtrain(cvp.training(k),:);
        y = ytrain(cvp.training(k),:);
        Xvalid = Xtrain(cvp.test(k),:);
        yvalid = ytrain(cvp.test(k),:);

        % Do the NCA for feature selectiom
        nca = fscnca(X,y,'FitMethod','exact', ...
             'Solver','sgd','Lambda',lambdavals(i), ...
             'IterationLimit',30,'GradientTolerance',1e-4, ...
             'Standardize',true);
                  
        % Calculate the loss values
        lossvals(i,k) = loss(nca,Xvalid,yvalid,'LossFunction','classiferror');
    end
end

% Average loss across folds for each lambda
meanloss = mean(lossvals,2);
[~,idx] = min(meanloss); % Find the index
bestlambda = lambdavals(idx); % best lambda

if 0
figure()
plot(lambdavals,meanloss,'ro-')
xlabel('Lambda')
ylabel('Loss (MSE)')
grid on
end

% fit nca on all training data
nca = fscnca(Xtrain,ytrain,'FitMethod','exact','Solver','sgd',...
    'Lambda',bestlambda,'Standardize',true,'Verbose',0);


if 0
figure()
plot(nca.FeatureWeights,'ro')
xlabel('Feature index')
ylabel('Feature weight')
grid on
xticks(1:length(features))
xticklabels(features)
end

tol    = 0.02; % tolerance in matlab example
selidx = find(nca.FeatureWeights > tol*max(nca.FeatureWeights)); % find features whose weights beat tolerance

% restrict features to this set
includedPredictorNames = predictorNames(selidx);
predictors = predictors(:,selidx);
isCategoricalPredictor = repmat(false,1,length(selidx));

%% PCA
%
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 'inf' values have to be treated as missing data for PCA.
numericPredictors(isinf(numericPredictors)) = NaN;

% Normalize the predictors prior to PCA
normalizationFcn = @(x) (x-nanmean(x,1))./nanstd(x,[],1);
normalizedPredictors = normalizationFcn(numericPredictors);
% Do PCA
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    normalizedPredictors);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = pca_perc/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep)), predictors(:, isCategoricalPredictor)];
%}

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
%predictFcn = @(x) oldPredictFcn(featureSelectionFcn(predictorExtractionFcn(x)));
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



