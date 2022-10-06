function [trueClass,predClass,AUC] = more_general_train_test(T,N,perc_train,model,predictors,response)

trueClass = cell(N,1);
predClass = cell(N,1);
AUC = nan(N,1);
npts = size(T,1);
all_pts = 1:npts;
ntrain = round(npts*perc_train);


for ib = 1:N
    while 1
        %% Do a random testing/training split
        r = randsample(npts,ntrain);
        training = ismember(all_pts,r);
        testing = ~ismember(all_pts,r);

        assert(isempty(intersect(find(training),find(testing))))

        Ttrain = T(training,:);
        Ttest = T(testing,:);

        %if length(unique(Ttest.(response))) == 1, continue; end % try again if only 1 class
        if length(unique(Ttest.(response))) == 1, break; end 

        %% Train the model
        tc = model(Ttrain,response,predictors);



        %% Predict the model on the testing data
        if isfield(tc,'predictFcn')
            predClass{ib} = tc.predictFcn(Ttest);
        else
            predClass{ib} = predict(tc,Ttest);
        end

        trueClass{ib} = Ttest.(response);

        if ~any(trueClass{ib}(~isnan(predClass{ib}))==1) || ~any(trueClass{ib}(~isnan(predClass{ib}))==0)
            continue; 
        end % try again if no true labels or no false labels
        [~,~,~,AUC(ib)] = perfcurve(trueClass{ib},predClass{ib},1);
        break
    end

end

end