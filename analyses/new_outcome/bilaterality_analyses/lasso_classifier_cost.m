function trainedClassifier = lasso_classifier_cost(trainingData,features,respVar,pc_perc,classes)

% this code trains a binary classifier. It first does PCA to reduce
% dimensionality, then does matlab's fitclinear with LASSO regularization
% and applies a cost function that is inversely weighted by the number of
% samples in the training data for each class.

% Seed a random number generator
rng(0)

%% Data prep
inputTable = trainingData;
predictorNames = features;
predictors = inputTable(:, predictorNames); % get the predictors
old_predictors = predictors;
response = inputTable.(respVar); % get the response


%% PCA
if length(features) ==1 
    % skip PCA
else
    numericPredictors = table2array(varfun(@double, predictors));
    numericPredictors(isinf(numericPredictors)) = NaN;
    
    % Do PCA with normalization
    w = 1./std(numericPredictors,[],1,"omitnan");
    [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
        numericPredictors,'centered',true,'VariableWeights',w);
    % Keep enough components to explain the desired amount of variance.
    explainedVarianceToKeepAsFraction = pc_perc/100;
    numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
    pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
    predictors = array2table(pcaScores(:,1:numComponentsToKeep));
end


%% Make a function to turn 1s and 0s into classes
response_bin = nan(length(response),1);
possible_responses = classes;
assert(length(possible_responses)==2) % confirm just 2
response_bin(strcmp(response,possible_responses{2})) = 1;
response_bin(strcmp(response,possible_responses{1})) = 0;
assert(sum(isnan(response_bin))==0)

% Make function to extract class back from one or zero
class_from_bin = @(x) possible_responses(x+1);


%% Establish cost function
n_in_class = nan(length(classes),1);
for i = 1:length(classes)
    n_in_class(i) = sum(strcmp(response,classes{i})); % get number in each class
end
cost = [0 n_in_class(2);n_in_class(1) 0];
%cost = [0 1;1 0];

%% Do the classifier
% default lambda for lasso is 1/n where n is training sample size, and this
% appears to be the lamda it is using in my dataset. No cross validation to
% find optimal lambda. I could potentially optimize this.

classifier = fitclinear(predictors,response,'Learner','logistic',...
    'ClassNames',classes,'cost',cost,'Regularization','Lasso');

% get the coefficients for the LR model
coef = [classifier.Bias; classifier.Beta];

% Create logit function with the trained coefficients (for the testing
% data)
lr_prob = @(x) glmval(coef,x,'logit');
lr_prediction = @(x) (x>=0.5); % make corresponding prediction rule with threshold 0.5

% alternate classifier
alt_class = fitcsvm(table2array(old_predictors),response,'KernelFunction','linear',...
    'ClassNames',classes,'cost',cost);


% Make final output anonumous functions
predictorExtractionFcn = @(t) t(:, predictorNames); % get the predictors

if length(features) == 1
    predictFcn = @(x) class_from_bin(lr_prediction(lr_prob(table2array(varfun(@double, predictorExtractionFcn(x)))))); % generate final class preciction
    
    probabilityFcn = @(x) lr_prob(table2array(varfun(@double, predictorExtractionFcn(x)))); % generate probability
    predFcn2 = @(x) classifier.predict(predictorExtractionFcn(x));
    guessScoreFcn = @(x) 1./(1+exp(-(table2array(varfun(@double, predictorExtractionFcn(x)))*coef(2:end)+coef(1))));
else
    pcaTransformationFcn = @(x) (table2array(varfun(@double, x)) - pcaCenters) .* w * pcaCoefficients; % apply PCA using the coefficients obtained from the training data
    invTransformationFcn = @(x) x'/pcaCoefficients/w+pcaCenters; % I only need this for understanding feature importance
    predictFcn = @(x) class_from_bin(lr_prediction(lr_prob(pcaTransformationFcn(predictorExtractionFcn(x))))); % generate final class preciction
    probabilityFcn = @(x) lr_prob(pcaTransformationFcn(predictorExtractionFcn(x))); % generate probability
    predFcn2 = @(x) classifier.predict(pcaTransformationFcn(predictorExtractionFcn(x)));
    guessScoreFcn = @(x) 1./(1+exp(-(pcaTransformationFcn(predictorExtractionFcn(x))*coef(2:end)+coef(1))));
end
altPredictFcn = @(x) alt_class.predict(table2array(predictorExtractionFcn(x)));
guessPred = @(x) guess_prediction_fcn(guessScoreFcn(x),classes);


% Add additional fields to the result struct
trainedClassifier.predictFcn = predictFcn;
trainedClassifier.RequiredVariables = predictorNames;
trainedClassifier.coef = coef;

if length(features) == 1
    trainedClassifier.PCACenters = [];
    trainedClassifier.PCACoefficients = [];
    trainedClassifier.pcaTransformationFcn = [];
    trainedClassifier.invTransformationFcn = [];
    trainedClassifier.pcaWeights = [];
else
    trainedClassifier.PCACenters = pcaCenters;
    trainedClassifier.PCACoefficients = pcaCoefficients;
    trainedClassifier.pcaTransformationFcn = pcaTransformationFcn;
    trainedClassifier.invTransformationFcn = invTransformationFcn;
    trainedClassifier.pcaWeights = w;
end

trainedClassifier.probabilityFcn = probabilityFcn;
trainedClassifier.predictorExtractionFcn = predictorExtractionFcn;

trainedClassifier.classes = classes;
trainedClassifier.n_in_class = n_in_class;
trainedClassifier.cost = cost;
trainedClassifier.altPredictFcn = altPredictFcn;
trainedClassifier.alt_class = alt_class;
trainedClassifier.classifier = classifier;
trainedClassifier.lr_prob = lr_prob;
trainedClassifier.coef = coef;
trainedClassifier.predFcn2 = predFcn2;
trainedClassifier.guessScoreFcn = guessScoreFcn;
trainedClassifier.guessPred = guessPred;



end

function pred= guess_prediction_fcn(x,classes)

% assume class 2 is positive and class 1 negative - double check
pos_class = classes{2};
neg_class = classes{1};

pred = cell(length(x),1);
for i = 1:length(x)

    if x(i) >= 0.5
        pred{i} = pos_class;
    else
        pred{i} = neg_class;
    end

end

end

