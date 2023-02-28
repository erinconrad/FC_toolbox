function out = classifier_wrapper(T,features,pca_perc,combine_br,just_spikes,rm_non_temporal)

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
npts = size(T,1);

% Remove non temporal patients if desired
if rm_non_temporal == 1
    temporal = strcmp(T.soz_locs,'temporal');
    T(~temporal,:) = [];
    npts = size(T,1);
elseif rm_non_temporal == 2 % only include non-temporal (excludes diffuse and multifocal)
    extra = strcmp(T.soz_locs,'other cortex') | strcmp(T.soz_locs,'frontal');
    T(~extra,:) = [];
    npts = size(T,1);
end

% Combine right and bilateral or left and bilateral
if combine_br == 1
    T.soz_lats(strcmp(T.soz_lats,'right') | strcmp(T.soz_lats,'bilateral')) = {'br'};
elseif combine_br == 2
    T.soz_lats(strcmp(T.soz_lats,'left') | strcmp(T.soz_lats,'bilateral')) = {'bl'};
end

% Initialize ROC parameters
classes = unique(T.(response));
all_scores = nan(npts,1);
nclasses = length(classes);
C = zeros(nclasses,nclasses); % left, right, bilateral
all_pred = cell(npts,1);

%% Do leave-one-patient-out classifier to predict laterality
for i = 1:npts
    i
    % split into training and testing
    Ttrain = T([1:i-1,i+1:end],:); % training all but current
    Ttest = T(i,:); % testing current

    % make sure they're distinct
    assert(isempty(intersect(Ttrain.names,Ttest.names)))

    % perform imputation of missing data
    for j = 1:size(Ttrain,2)
        a = Ttrain{:,j};
        if ~isnumeric(a), continue; end

        a(isnan(a)|abs(a)<1e-10) = nanmedian(a);
        Ttrain{:,j} = a;

        b = Ttest{:,j};
        b(isnan(b)|abs(b)<1e-10) = nanmedian(a); % impute with training data median
        Ttest{:,j} = b;
    end

    % Dumb spikes - binarize spikes according to which side has more
    if just_spikes == 2
        for j = 1:size(Ttrain,2)
            if ~isnumeric(Ttrain{:,j}), continue; end
            Ttrain{Ttrain{:,j}>0,j} = 1; Ttrain{Ttrain{:,j}<0,j} = -1;
            Ttest{Ttest{:,j}>0,j} = 1; Ttest{Ttest{:,j}<0,j} = -1;
        end

    end

    % Do the lasso calculator
    tc = lasso_classifier(Ttrain,features,response,pca_perc,classes);

    % Get score
    all_scores(i) = tc.probabilityFcn(Ttest);

    % make prediction on left out
    pred = tc.predictFcn(Ttest);
    all_pred{i} = pred{1};
    
    % compare to true
    true = Ttest.(response);

    % which row to add to confusion matrix (the true value)
    which_row = find(strcmp(true,classes));

    % which column to add to confusion matrix (the predicted value)
    which_column = find(strcmp(pred,classes));

    C(which_row,which_column) = C(which_row,which_column) + 1;


end

% Prepare output structure
out.scores = all_scores;
out.class = T.(response);
out.pos_class = classes{2};
out.all_pred = all_pred;
out.C = C;
out.unique_classes = classes;
out.npts = npts;

end