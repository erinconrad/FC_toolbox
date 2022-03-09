function soz_roc_out = classify_soz(for_bootstrap)

%{
Logistic regression model. I confirmed that when I take the training and
testing data and put them into Matlab's classification learner with an LR
model it gives the same AUC

%}

%% Parameters
include_post  = 0;
prop_train = 2/3;
do_norm = 1;

colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];

%% Seed rng (for splitting testing and training data)
% If I don't seed this the AUC usually bounces around 0.73-0.82
if for_bootstrap == 0 % if 0, I seed the rng to get the ORs. If 1, I don't seed to get stats
    rng(0)
end

locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);


%% Get stuff
rate_sw = out.bin_out.all_elecs_rates_sw; %rates wake and sleep
%rl_sw = out.bin_out.all_elecs_rl_sw;
rate_all = out.bin_out.all_elecs_rates;
%rl_all = out.bin_out.all_elecs_rl;
soz = out.bin_out.all_is_soz;
pt_idx = (1:length(rate_sw))';
rate_post = out.sz_out.post_ictal_rates;

%% Put data into friendly format
vec_rate_sleep = [];
vec_rate_wake = [];
vec_rate_all =[];
vec_rate_post = [];
%vec_rl_sleep = [];
%vec_rl_wake = [];
%vec_rl_all = [];
vec_pt_idx = [];
vec_soz = [];
for ip = 1:length(pt_idx)
    
    curr_rate_sw = rate_sw{ip};
    curr_rate_all = rate_all{ip};
    %curr_rl_sw = rl_sw{ip};
    %curr_rl_all = rl_all{ip};
    curr_rate_post = rate_post{ip};
    
    if do_norm
        % normalize the rate across electrodes (this is so that patients
        % with higher spike rates in general will get the same weight as
        % patients with lower spike rates)
        curr_rate_sw = (curr_rate_sw - nanmean(curr_rate_sw,1))./...
            nanstd(curr_rate_sw,[],1);
        curr_rate_all = (curr_rate_all-nanmean(curr_rate_all))./...
            nanstd(curr_rate_all);
        %{
        curr_rl_sw = (curr_rl_sw - nanmean(curr_rl_sw,1))./...
            nanstd(curr_rl_sw,[],1);
        curr_rl_all = (curr_rl_all-nanmean(curr_rl_all))./...
            nanstd(curr_rl_all);
            %}
        
        curr_rate_post = (curr_rate_post - nanmean(curr_rate_post))./...
            nanstd(curr_rate_post);
    end
    
    vec_rate_sleep = [vec_rate_sleep;curr_rate_sw(:,2)];
    vec_rate_wake = [vec_rate_wake;curr_rate_sw(:,1)];
    vec_rate_all = [vec_rate_all;curr_rate_all];
    vec_rate_post = [vec_rate_post;curr_rate_post];
    
    %vec_rl_sleep = [vec_rl_sleep;curr_rl_sw(:,2)];
    %vec_rl_wake = [vec_rl_wake;curr_rl_sw(:,1)];
    %vec_rl_all = [vec_rl_all;curr_rl_all];
    
    vec_pt_idx = [vec_pt_idx;repmat(pt_idx(ip),length(rate_sw{ip}(:,2)),1)];
    
    vec_soz = [vec_soz;soz{ip}'];
    
end



%% Make table
%{
T = table(vec_soz,vec_pt_idx,vec_rate_sleep,vec_rate_wake,...
    vec_rl_sleep,vec_rl_wake,vec_rate_all,vec_rl_all,vec_rate_post);
%}
T = table(vec_soz,vec_pt_idx,vec_rate_sleep,vec_rate_wake);
T.vec_pt_idx = (T.vec_pt_idx);

%% Remove nan and inf rows
nan_rows = isnan(T.vec_rate_sleep) | isnan(T.vec_rate_wake) | isnan(T.vec_soz) | isnan(T.vec_pt_idx);
T(nan_rows,:) = [];
%vec_pt_idx(nan_rows) = [];


%% Divide into training and testing set
training = randsample(length(pt_idx),floor(length(pt_idx)*prop_train));
training_idx = ismember(T.vec_pt_idx,training);
testing_idx = ~ismember(T.vec_pt_idx,training);

T_train = T(training_idx,:);
T_test = T(testing_idx,:);

% Confirm that I separated by patients
train_pts = unique(T_train.vec_pt_idx);
test_pts = unique(T_test.vec_pt_idx);
assert(isempty(intersect(train_pts,test_pts)))
%assert(isequal(str2double(string(sort([train_pts;test_pts]))),(1:96)'))
% this above statement is not true because pt 7 has only one sleep period
% (with all nans for spike rate, maybe bad time...)
assert(sum(isnan(table2array(T)),'all')==0)

%% Train model with glm model

glm = fitglm(T_train,...
    'vec_soz ~ vec_rate_sleep + vec_rate_wake',...
    'Distribution','Poisson','Link','log');

%{
if include_post
    glme = fitglme(T_train,...
        'vec_soz ~ vec_rate_sleep + vec_rate_wake + vec_rate_post + (1|vec_pt_idx)',...
        'Distribution','Poisson','Link','log');
else
    glme = fitglme(T_train,...
        'vec_soz ~ vec_rate_sleep + vec_rate_wake + (1|vec_pt_idx)',...
        'Distribution','Poisson','Link','log');
end
%}
%% For generating odds ratios, do a glme
T_train_glme = T_train;
T_train_glme.vec_pt_idx = nominal(T_train_glme.vec_pt_idx);
glme = fitglme(T_train,...
    'vec_soz ~ vec_rate_sleep + vec_rate_wake + (1|vec_pt_idx)',...
    'Distribution','Poisson','Link','log');



%% Test model on testing data
params = glm.CoefficientNames;
asum = zeros(size(T_test,1),1);
asum = asum + glm.Coefficients.Estimate(1);
for p = 2:length(params)
    est = glm.Coefficients.Estimate(p);
    asum = asum +  T_test.(params{p})*est;
end


classification = logistic(asum);
nll = NLL(T_test.vec_soz,classification);
pred_nll = log(2)*length(classification);
%fprintf('\nTest data log likelihood: %1.1f\npredicted by chance: %1.1f\n',nll,pred_nll);

% Calculate odds ratios and CI for odds ratios of each predictor
sleep_beta = glme.Coefficients{2,2};
wake_beta = glme.Coefficients{3,2};
sleep_se = glme.Coefficients{2,3};
wake_se = glme.Coefficients{3,3};
sleep_t = glme.Coefficients{2,4};
wake_t = glme.Coefficients{3,4};
sleep_p = glme.Coefficients{2,6};
wake_p = glme.Coefficients{3,6};


sleep_or = exp(sleep_beta);
wake_or = exp(wake_beta);
sleep_ci95 = [exp(sleep_beta - 1.96*sleep_se),exp(sleep_beta + 1.96*sleep_se)];
wake_ci95 = [exp(wake_beta - 1.96*wake_se),exp(wake_beta + 1.96*wake_se)];

all_soz = T_test.vec_soz==1;
all_no_soz = T_test.vec_soz==0;
class_soz = classification(all_soz);
class_no_soz = classification(all_no_soz);
[roc,auc,disc,disc_I] = calculate_roc(class_no_soz,class_soz,1e3);
soz_roc_out.roc = roc;
soz_roc_out.auc = auc;
soz_roc_out.glme = glme;
soz_roc_out.class_soz = class_soz;
soz_roc_out.class_no_soz = class_no_soz;
soz_roc_out.T_test = T_test;
soz_roc_out.T_train = T_train;
soz_roc_out.disc = disc;
soz_roc_out.disc_I = disc_I;
soz_roc_out.sleep_or = sleep_or;
soz_roc_out.wake_or = wake_or;
soz_roc_out.sleep_ci95 = sleep_ci95;
soz_roc_out.wake_ci95 = wake_ci95;
soz_roc_out.sleep_p = sleep_p;
soz_roc_out.wake_p = wake_p;
soz_roc_out.sleep_t = sleep_t;
soz_roc_out.wake_t = wake_t;


%{
%% Are patients with larger number of high probability classifications more likely to be multifocal or diffuse?
asum = zeros(size(T,1),1);
asum = asum + glme.Coefficients.Estimate(1);
for p = 2:length(params)
    est = glme.Coefficients.Estimate(p);
    asum = asum +  T.(params{p})*est;
end
classification = logistic(asum);

% Loop over patients

pts = unique(T.vec_pt_idx);
avg_class = nan(length(pts),1);
for i = 1:length(pts)
    p = pts(i);
    rows = find(T.vec_pt_idx==p);
    avg_class(i) = sum(classification(rows) > 0.1);
    %avg_class(i) = nanmean(classification(rows));
end
locs = circ_out.all_locs;
locs = locs(pts);
diffuse = strcmp(locs,'multifocal') | strcmp(locs,'diffuse');
if 1
    figure
    plot(1+0.05*randn(sum(diffuse),1),avg_class(diffuse),'o')
    hold on
    plot(2+0.05*randn(sum(~diffuse),1),avg_class(~diffuse),'o')
end
%}

%% Alt ROC method
%
labels = cell(length(classification),1);
labels(all_soz) = {'SOZ'};
labels(all_no_soz) = {'Not SOZ'};
posclass = 'SOZ';
scores = classification;
[X,Y,T,AUC,opt] = perfcurve(labels,scores,posclass);
soz_roc_out.alt_auc = AUC;
soz_roc_out.X = X;
soz_roc_out.Y = Y;
%}

%% ROC
if 0
figure
plot(roc(:,1),roc(:,2),'k-','linewidth',2)
hold on
plot([0 1],[0 1],'k--','linewidth',2)
%plot(roc(disc_I,1),roc(disc_I,2),'*','markersize',15,'linewidth',2,'color',colors(5,:));
%text(roc(disc_I,1)+0.01,roc(disc_I,2)-0.05,'SOZ cutoff','fontsize',15,'color',colors(5,:));
xlabel('False positive rate')
ylabel('True positive rate')
legend(sprintf('AUC %1.2f',auc),'location','southeast','fontsize',15)
set(gca,'fontsize',15)
title('SOZ identification accuracy')
end

end


function out = logistic(x)

out = 1./(1+exp(-x));

end

function nll = NLL(X,P)

nll = nansum(-X.*log(P) - (1-X).*log(1-P));
end