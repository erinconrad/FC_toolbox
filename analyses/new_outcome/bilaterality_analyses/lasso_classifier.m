function trainedClassifier = lasso_classifier(trainingData,features,respVar,pc_perc,classes)

%% Data prep
inputTable = trainingData;
predictorNames = features;
predictors = inputTable(:, predictorNames);
response = inputTable.(respVar);


%% PCA
numericPredictors = table2array(varfun(@double, predictors));
numericPredictors(isinf(numericPredictors)) = NaN;

% Do PCA
w = 1./std(numericPredictors,[],1,"omitnan");
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors,'centered',true,'VariableWeights',w);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = pc_perc/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = array2table(pcaScores(:,1:numComponentsToKeep));

%% LASSO logistic regression
% Make response ones and zeros
response_bin = nan(length(response),1);
possible_responses = classes;
assert(length(possible_responses)==2) % confirm just 2
response_bin(strcmp(response,possible_responses{2})) = 1;
response_bin(strcmp(response,possible_responses{1})) = 0;
assert(sum(isnan(response_bin))==0)

% Make function to extract class back from one or zero
class_from_bin = @(x) possible_responses(x+1); % if 0->possible_responses{1}, if 1->possible_responses{2}
assert(isequal(class_from_bin(response_bin),response)) % confirm I get original classes back

% do the Lasso
[B,FitInfo] = lassoglm(table2array(predictors),response_bin,'binomial','CV',5);
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
B0 = FitInfo.Intercept(idxLambdaMinDeviance);
coef = [B0; B(:,idxLambdaMinDeviance)];

% Do logistic regression
lr_prob = @(x) glmval(coef,x,'logit');
lr_prediction = @(x) (x>=0.5);


%% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
pcaTransformationFcn = @(x) (table2array(varfun(@double, x)) - pcaCenters) .* w * pcaCoefficients;
predictFcn = @(x) class_from_bin(lr_prediction(lr_prob(pcaTransformationFcn(predictorExtractionFcn(x)))));
probabilityFcn = @(x) lr_prob(pcaTransformationFcn(predictorExtractionFcn(x)));

% Add additional fields to the result struct
trainedClassifier.predictFcn = predictFcn;
trainedClassifier.RequiredVariables = predictorNames;
trainedClassifier.PCACenters = pcaCenters;
trainedClassifier.coef = coef;
trainedClassifier.PCACoefficients = pcaCoefficients;
trainedClassifier.probabilityFcn = probabilityFcn;
trainedClassifier.pcaTransformationFcn = pcaTransformationFcn;
trainedClassifier.predictorExtractionFcn = predictorExtractionFcn;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2022a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');



