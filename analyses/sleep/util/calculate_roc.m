function [roc,auc,disc] = calculate_roc(array1,array2,npts)

% positive = I say it is sleep
% negative = I say it is wake
% "disease" = really is asleep
% FP = I say sleep but really awake
% TN = I say wake, really wake
% FN = I say wake, really asleep
% TP = I say sleep, really asleep

%% Parameters

%% Initialize output
roc = nan(npts,2); % 1 = FPR, 2 = TPR

%% Get the range of potential values
min_val = min([array1;array2]);
max_val = max([array1;array2]);

%% Create steps for values
test_vals = linspace(min_val,max_val,npts);

%% Loop over test vals and calculate stats
for i = 1:npts
    val = test_vals(i); % the threshold a/d below which I am calling it sleep
    
    %% Calculate FPR (FP/(FP+TN))
    FP = sum(array2 < val); % number of waking cases with low alpha/delta (and so I call sleep)
    TN = sum(array2 >= val); % number of waking cases with high alpha/delta
    FPR = FP/(FP+TN);
    
    %% Calculate TPR (TP/(TP+FN))
    TP = sum(array1 < val); % number of sleeping cases with low alpha delta
    FN = sum(array1 >= val); % number of sleeping cases with high alpha delta
    TPR = TP/(TP+FN);
    
    roc(i,:) = [FPR TPR];
    
end

%% Calculate AUC
auc = trapz(roc(:,1),roc(:,2));

%% Find discriminant
% Find the threshold normalized alpha delta ratio that best separates wake
% from sleep
goodness = (ones(size(roc,1),1) - roc(:,1)) + roc(:,2); % 1-FPR + TPR;
[~,disc] = max(goodness);
disc = test_vals(disc);

end