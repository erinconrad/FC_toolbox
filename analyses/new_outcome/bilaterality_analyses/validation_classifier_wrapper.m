function out = validation_classifier_wrapper(T,train,test,features,pca_perc,...
    combine_br,just_spikes,rm_non_temporal,which_ref)
% this code is the wrapping function to train a model on a training dataset
% and test on an external validation set

% Define response
response = 'soz_lats';

% Restrict to spike features if desired
spike_features = features(contains(features,'spikes') & contains(features,which_ref));

if just_spikes == 1 || just_spikes == 2
    features = spike_features;
end

% Remove patients without response
empty_class = cellfun(@isempty,T.(response));
assert(sum(empty_class)==0) % I should already have removed these
%T(empty_class,:) = [];
%train(empty_class) = [];
%test(empty_class) = [];
npts = size(T,1);

% Remove non temporal patients if desired
if rm_non_temporal == 1
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];
    train(~temporal) = []; test(~temporal) = [];
    npts = size(T,1);
elseif rm_non_temporal == 2 % only include non-temporal (excludes diffuse and multifocal)
    extra = strcmp(T.soz_locs,'other cortex') | strcmp(T.soz_locs,'frontal');
    T(~extra,:) = [];
    train(~extra) = []; test(~extra) = [];
    npts = size(T,1);
end

% If I am doing a single feature model, exclude patient(s) with nan for
% variable
if just_spikes == 1 || just_spikes == 2   
    nan_feature = isnan(T{:,spike_features});
    out.rm_nan_spikes = sum(nan_feature);
    T(nan_feature,:) = [];
    train(nan_feature) = []; test(nan_feature) = [];
    npts = size(T,1);
end



% Combine right and bilateral or left and bilateral
if combine_br == 1
    T.soz_lats(strcmp(T.soz_lats,'right') | strcmp(T.soz_lats,'bilateral')) = {'br'};
elseif combine_br == 2
    T.soz_lats(strcmp(T.soz_lats,'left') | strcmp(T.soz_lats,'bilateral')) = {'bl'};
elseif combine_br == 3
    T.soz_lats(strcmp(T.soz_lats,'left') | strcmp(T.soz_lats,'right')) = {'lr'};
end

% Initialize ROC parameters
classes = unique(T.(response));
nclasses = length(classes);
C = zeros(nclasses,nclasses); % left, right, bilateral

% Establish training and testing data
Ttrain = T(train,:);
Ttest = T(test,:);

% Perform imputation of missing data
for j = 1:size(Ttrain,2)
    a = Ttrain{:,j};
    if ~isnumeric(a), continue; end

    a(isnan(a)) = nanmedian(a);
    Ttrain{:,j} = a;

    b = Ttest{:,j};
    b(isnan(b)) = nanmedian(a); % impute with training data median
    Ttest{:,j} = b;
end

% Binarize spikes according to which side has more
if just_spikes == 2
    binTtrain = Ttrain;
    binTtest = Ttest;

    for j = 1:size(Ttrain,2)
        if ~isnumeric(Ttrain{:,j}), continue; end
        binTtrain{Ttrain{:,j}>0,j} = 1; binTtrain{Ttrain{:,j}<0,j} = 0; % 1 if positive, 0 otherwise
        binTtest{Ttest{:,j}>0,j} = 1; binTtest{Ttest{:,j}<0,j} = 0;
        binTtrain{:,j} = logical(binTtrain{:,j});
        binTtest{:,j} = logical(binTtest{:,j});
    end
    Ttrain = binTtrain;
    Ttest = binTtest;
end

if combine_br == 0
    tc = general_classifier(Ttrain,'svm_cost',features,response,pca_perc,1e2,length(features),classes);
    all_scores = [];
else
    % train the model on the training data
    tc = lasso_classifier_cost(Ttrain,features,response,pca_perc,classes);

    % get the probability of the testing data
    all_scores = tc.probabilityFcn(Ttest);
end
all_names = Ttest.names;
%all_pred = tc.predictFcn(Ttest);
%all_pred = tc.classifier.predict(Ttest);
%alt_pred = tc.altPredictFcn(Ttest);
all_pred = tc.predFcn2(Ttest);
alt_pred = tc.guessPred(Ttest);

if contains(features,'samp')
    % for some reason, the class labels get flipped in the subsampling
    % analysis so don't make this requirement

else
    assert(isequal(all_pred,alt_pred))
end

%% Re-derive feature weights
coef = tc.coef(2:end);
pcaCoefficients = tc.PCACoefficients;
pcaCenters = tc.PCACenters;
w = tc.pcaWeights;
if length(coef) > 1
    coef_feature_space = tc.invTransformationFcn(coef);
    [~,I] = sort(abs(coef_feature_space),'descend');
    sorted_features = features(I);
    out.sorted_features = sorted_features;
    out.coefs = coef_feature_space(I);
end

%% Make confusion matrix
response_true = Ttest.(response);
%{
positive = strcmp(Ttest.(response),classes{2});
negative = strcmp(Ttest.(response),classes{1});
pred_positive = strcmp(all_pred,classes{2});
pred_negative = strcmp(all_pred,classes{1});
C(1,1) = sum(positive & pred_positive);
C(1,2) = sum(positive & pred_negative);
C(2,1) = sum(negative & pred_positive);
C(2,2) = sum(negative & pred_negative);
%}
for i = 1:length(all_pred)
    pred = all_pred{i};
    actual = response_true{i};

    % which row to add to confusion matrix (the true value)
    which_row = find(strcmp(actual,classes));

    % which column to add to confusion matrix (the predicted value)
    which_column = find(strcmp(pred,classes));

    C(which_row,which_column) = C(which_row,which_column) + 1;

end

%% Prepare output structure
out.scores = all_scores;
out.class = Ttest.(response);
out.pos_class = classes{2};
out.all_pred = all_pred;
out.C = C;
out.unique_classes = classes;
out.npts = npts;
out.names = all_names;
out.tc = tc;
%out.alt_pred = alt_pred;

out.accuracy = (C(1,1)+C(2,2))/sum(C,'all');

recall = nan(nclasses,1);
for i = 1:nclasses
    tp = C(i,i);
    fn = sum(C(i,~ismember(1:nclasses,i))); 
    recall(i) = tp/(tp+fn); % tp is number correctly predicted to be in class, tp + fn is everyone in the class
end
out.balanced_accuracy = mean(recall);


end