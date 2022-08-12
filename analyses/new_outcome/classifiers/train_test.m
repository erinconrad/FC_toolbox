function [confusion,test_accuracy,cats] = train_test(T,N,perc_train,model,forn)

npts = size(T,1);
if isempty(N)
    % do LOO
    loo = 1;
    N = npts;
else
    loo = 0;
end

%% Initialize confusion matrix
cats = unique(T.SOZ);
ncats = length(cats);
confusion = zeros(ncats,ncats);

%confusion = cell(N,1);
test_accuracy = nan(N,1);

all_pts = 1:npts;
ntrain = round(npts*perc_train);

for ib = 1:N
    %% Do a random testing/training split
    if loo == 1
        % LOO
        testing = ismember(all_pts,ib);
        training = ~ismember(all_pts,ib);
    else
        r = randsample(npts,ntrain);
        training = ismember(all_pts,r);
        testing = ~ismember(all_pts,r);
    end

    assert(isempty(intersect(find(training),find(testing))))

    Ttrain = T(training,:);
    Ttest = T(testing,:);

    %% Train the model
    tc = model(Ttrain,forn);

    %% Predict the model on the testing data
    if isfield(tc,'predictFcn')
        predClass = tc.predictFcn(Ttest);
    else
        predClass = predict(tc,Ttest);
    end

    if isequal(model,@outcome_logistic_regression)
        predClassOld = predClass;
        predClass = cell(length(predClassOld),1);
        predClass(predClassOld > 0.5) = {'good'};
        predClass(predClassOld <= 0.5) = {'bad'};
        trueClassOld = Ttest.outcome;
        trueClass = cell(length(trueClassOld),1);
        trueClass(trueClassOld > 0.5) = {'good'};
        trueClass(trueClassOld <= 0.5) = {'bad'};

    elseif isequal(model,@sozTree) || isequal(model,@sozTreePCA)
        trueClass = Ttest.SOZ;
    end


    %% add things to appropriate parts of confusion matrix
    for ic = 1:ncats
        for jc = 1:ncats
            confusion(ic,jc) = confusion(ic,jc) + ...
                sum(strcmp(trueClass,cats{ic}) & strcmp(predClass,cats{jc}));
        end
    end

    test_accuracy(ib) = sum(cellfun(@(x,y) strcmp(x,y),trueClass,predClass))./length(predClass);

    
end



end