function soz_roc_out = classify_soz

%{
mixed effects Logistic regression model
Need to check this!!!
%}

%% Parameters
include_sw = 1;
include_post  = 1;
prop_train = 2/3;
do_norm = 1;

colors = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.4660, 0.6740, 0.1880;...
    0.4940, 0.1840, 0.5560;...
    0.6350, 0.0780, 0.1840];

%% Seed rng (because it's a MC test)
rng(0)

locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

%% Load out file and get roc stuff
out = load([out_folder,'out.mat']);
out = out.out;

%% Unpack substructures
unpack_any_struct(out);
out_folder = [results_folder,'analysis/sleep/'];

%% Get stuff
rate_sw = out.bin_out.all_elecs_rates_sw; %rates wake and sleep
rl_sw = out.bin_out.all_elecs_rl_sw;
rate_all = out.bin_out.all_elec_rates;
rl_all = out.bin_out.all_elecs_rl;
soz = out.bin_out.all_is_soz;
pt_idx = (1:length(rate_sw))';
rate_post = out.sz_out.post_ictal_rates;

%% Put data into friendly format
vec_rate_sleep = [];
vec_rate_wake = [];
vec_rate_all =[];
vec_rate_post = [];
vec_rl_sleep = [];
vec_rl_wake = [];
vec_rl_all = [];
vec_pt_idx = [];
vec_soz = [];
for ip = 1:length(pt_idx)
    
    curr_rate_sw = rate_sw{ip};
    curr_rate_all = rate_all{ip};
    curr_rl_sw = rl_sw{ip};
    curr_rl_all = rl_all{ip};
    curr_rate_post = rate_post{ip};
    
    if do_norm
        curr_rate_sw = (curr_rate_sw - nanmedian(curr_rate_sw,1))./...
            iqr(curr_rate_sw,1);
        curr_rate_all = (curr_rate_all-nanmedian(curr_rate_all))./...
            iqr(curr_rate_all);
        
        curr_rl_sw = (curr_rl_sw - nanmedian(curr_rl_sw,1))./...
            iqr(curr_rl_sw,1);
        curr_rl_all = (curr_rl_all-nanmedian(curr_rl_all))./...
            iqr(curr_rl_all);
        
        curr_rate_post = (curr_rate_post - nanmedian(curr_rate_post))./...
            iqr(curr_rate_post);
    end
    
    vec_rate_sleep = [vec_rate_sleep;curr_rate_sw(:,2)];
    vec_rate_wake = [vec_rate_wake;curr_rate_sw(:,1)];
    vec_rate_all = [vec_rate_all;curr_rate_all];
    vec_rate_post = [vec_rate_post;curr_rate_post];
    
    vec_rl_sleep = [vec_rl_sleep;curr_rl_sw(:,2)];
    vec_rl_wake = [vec_rl_wake;curr_rl_sw(:,1)];
    vec_rl_all = [vec_rl_all;curr_rl_all];
    
    vec_pt_idx = [vec_pt_idx;repmat(pt_idx(ip),length(rl_sw{ip}(:,2)),1)];
    
    vec_soz = [vec_soz;soz{ip}'];
    
end



%% Make table
T = table(vec_soz,vec_pt_idx,vec_rate_sleep,vec_rate_wake,...
    vec_rl_sleep,vec_rl_wake,vec_rate_all,vec_rl_all,vec_rate_post);
T.vec_pt_idx = nominal(T.vec_pt_idx);

%% Divide into training and testing set
training = randsample(length(pt_idx),floor(length(pt_idx)*prop_train));
training_idx = ismember(vec_pt_idx,training);
testing_idx = ~ismember(vec_pt_idx,training);

T_train = T(training_idx,:);
T_test = T(testing_idx,:);

%% Train model with glme model
if include_post
    glme = fitglme(T_train,...
        'vec_soz ~ vec_rate_sleep + vec_rate_wake + vec_rate_post + (1|vec_pt_idx)',...
        'Distribution','Poisson','Link','log');
else
    glme = fitglme(T_train,...
        'vec_soz ~ vec_rate_sleep + vec_rate_wake + (1|vec_pt_idx)',...
        'Distribution','Poisson','Link','log');
end



%% Test model on testing data
params = glme.CoefficientNames;
asum = zeros(size(T_test,1),1);
asum = asum + glme.Coefficients.Estimate(1);
for p = 2:length(params)
    est = glme.Coefficients.Estimate(p);
    asum = asum +  T_test.(params{p})*est;
end


classification = logistic(asum);
nll = NLL(T_test.vec_soz,classification);
pred_nll = log(2)*length(classification);
fprintf('\nTest data log likelihood: %1.1f\npredicted by chance: %1.1f\n',nll,pred_nll);

all_soz = T_test.vec_soz==1;
all_no_soz = T_test.vec_soz==0;
class_soz = classification(all_soz);
class_no_soz = classification(all_no_soz);
[roc,auc,disc,disc_I] = calculate_roc(class_no_soz,class_soz,1e3);
soz_roc_out.roc = roc;
soz_roc_out.auc = auc;
soz_roc_out.glme = glme;

%% ROC
if 1
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