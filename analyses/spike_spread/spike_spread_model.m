function spike_spread_model

%% Seed rng
rng(0);
prop_train = 2/3;

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
int_folder = [results_folder,'analysis/intermediate/'];
out_folder = [results_folder,'analysis/spike_spread/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Load out file
out = load([out_folder,'out.mat']);
out = out.out;

%% table
T = out.T;

%% randomly divide patients into testing and training
npts = length(unique(T.all_pt_idx));
training = randsample(npts,floor(npts*prop_train));
training_idx = ismember(T.all_pt_idx,training);
testing_idx = ~ismember(T.all_pt_idx,training);

T_train = T(training_idx,:);
T_test = T(testing_idx,:);

% Confirm that I separated by patients
train_pts = unique(T_train.all_pt_idx);
test_pts = unique(T_test.all_pt_idx);
assert(isempty(intersect(train_pts,test_pts)))
assert(isequal(str2double(string(sort([train_pts;test_pts]))),(1:96)'))

%% Remove nan and inf rows
nan_train = isnan(T_train.all_coaw) | isnan(T_train.all_fcw) | isnan(T_train.all_distw);
T_train(nan_train,:) = [];

%% Train model
lm = fitlm(T,'all_coaw ~ all_fcw + all_distw');

%% Test model on testing data
params = lm.CoefficientNames;
asum = zeros(size(T_test,1),1);
asum = asum + lm.Coefficients.Estimate(1);
for p = 2:length(params)
    est = lm.Coefficients.Estimate(p);
    asum = asum +  T_test.(params{p})*est;
end

resid = (T_test.all_coaw - asum);
sum_sq_res = nansum(resid.^2);

%% Compare simpler model just distance
lm_simple = fitlm(T,'all_coaw ~ all_distw');

%% Test model on testing data
params = lm_simple.CoefficientNames;
asum = zeros(size(T_test,1),1);
asum = asum + lm_simple.Coefficients.Estimate(1);
for p = 2:length(params)
    est = lm_simple.Coefficients.Estimate(p);
    asum = asum +  T_test.(params{p})*est;
end

resid = (T_test.all_coaw - asum);
sum_sq_res_simple = nansum(resid.^2);

end