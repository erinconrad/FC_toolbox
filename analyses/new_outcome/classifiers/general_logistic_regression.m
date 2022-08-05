function all_auc = general_logistic_regression(T,params)

%% Remove nan rows
T(any(ismissing(T),2),:) = [];

loo = params.loo;
N = params.Nsplits;
perc_train = params.perc_train;
predictorNames = params.predictors;
response = params.response;
pt_id = T.(params.pt_id);
npts = length(unique(pt_id));
ntrain = round(perc_train*npts);

if loo == 1
    N = npts;
end
all_auc = nan(N,1);



for i = 1:N
    
    %% Do a random testing/training split
    if loo
        training = ~ismember(pt_id,i);
        testing = ismember(pt_id,i);
    else
        r = randsample(unique(pt_id),ntrain);
        training = ismember(pt_id,r);
        testing = ~ismember(pt_id,r);
    end

    assert(isempty(intersect(find(training),find(testing))))

    Ttrain = T(training,:);
    Ttest = T(testing,:);
    
    if isempty(Ttest) || sum(Ttest.SOZ) == 0
        continue
    end

    if ~loo
    assert(isequal(unique(Ttrain.(params.pt_id)),unique(r)))
    assert(isempty(intersect(Ttest.(params.pt_id),r)))
    end

    glm = fitglm(Ttrain,'predictorvars',predictorNames,'responsevar',response,'distribution','binomial');
    scores = predict(glm,Ttest);
    [~,~,~,AUC] = perfcurve(Ttest.SOZ,scores,1);
    all_auc(i) = AUC;
    
    
end

end