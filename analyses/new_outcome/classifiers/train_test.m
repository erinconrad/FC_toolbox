function [confusion,test_accuracy] = train_test(T,N,perc_train,model)

confusion = cell(N,1);
test_accuracy = nan(N,1);
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

        %% Train the model
        tc = model(Ttrain);

        %% Predict the model on the testing data
        predClass = tc.predictFcn(Ttest);
        trueClass = Ttest.SOZ;


        C = confusionmat(trueClass,predClass);
        accuracy = sum(cellfun(@(x,y) strcmp(x,y),trueClass,predClass))/numel(predClass);
        
        % Reject if funny errors
        if size(C,1) == 4
            continue
        else
            confusion{ib} = C;
            test_accuracy(ib) = accuracy;
            break
        end
            
    end
        
    
end

confusion = cat(3,confusion{:});


end