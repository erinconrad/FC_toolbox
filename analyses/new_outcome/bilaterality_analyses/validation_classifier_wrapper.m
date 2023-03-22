function out = validation_classifier_wrapper(T,train,test,features,pca_perc,combine_br,just_spikes,rm_non_temporal)

% Define response
response = 'soz_lats';

% Restrict to spike features if desired
spike_features = features(cellfun(@(x) contains(x,'spikes'),features));
if just_spikes == 1 || just_spikes == 2
    features = spike_features;
end

% Remove patients without response
empty_class = cellfun(@isempty,T.(response));
T(empty_class,:) = [];
train(empty_class) = [];
test(empty_class) = [];
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

%% Train the model on the training data
Ttrain = T(train,:);
Ttest = T(test,:);
tc = lasso_classifier(Ttrain,features,response,pca_perc,classes);
all_scores = tc.probabilityFcn(Ttest);
all_names = Ttest.names;
all_pred = tc.predictFcn(Ttest);


% Prepare output structure
out.scores = all_scores;
out.class = Ttest.(response);
out.pos_class = classes{2};
out.all_pred = all_pred;
out.C = C;
out.unique_classes = classes;
out.npts = npts;
out.names = all_names;
out.tc = tc;

end