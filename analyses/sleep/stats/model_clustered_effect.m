function [mean_estimate,ste_estimate] = model_clustered_effect(things_mat)

%% Reshape to a single vector
things = things_mat;
npts = size(things,1); % assumes patient is first dimenesion
pt_idx = repmat((1:npts)',1,size(things,2));

% vectorize
things = things(:);
pt_idx = pt_idx(:);

T = table(things,pt_idx);
T.pt_idx = nominal(T.pt_idx);
mdl = fitlme(T,'things~1+(1|pt_idx)');

%% Get intercept estimate and standard error
mean_estimate = mdl.Coefficients{1,2}; % intercept estimate
ste_estimate = mdl.Coefficients{1,3}; % intercept SE

% make sure the estimate is close to the unweighted mean
assert(abs(mean_estimate - nanmean(things)) < 5e-2);



end