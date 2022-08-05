function soz_predictor

%{
NEED TO PREP LOO AND COMPARE ACROSS OUTCOMES
%}

which_outcome = 'ilae';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;

 %% get variables of interest
ilae = data.all_two_year_ilae;
engel = data.all_two_year_engel;
surgery = data.all_surgery;
good_spikes = data.good_spikes;
rate = data.all_spikes;
soz = data.all_soz_bin;
locs = data.all_locs;
spike_rl = data.all_rl;
pc = data.all_fc;
coh = data.all_coh;
npts = length(rate);

%% Get outcome
switch which_outcome
    case 'ilae'
        outcome = ilae;
    case 'engel'
        outcome = engel;
end

%% Find good and bad outcome
outcome_num = cellfun(@(x) parse_outcome(x,which_outcome),outcome);

%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);

%% Find those with non-empty outcomes
non_empty = cellfun(@(x) ~isempty(x), outcome);

%% Find those with non empty outcomes, resection or ablation, and good spikes
complete = non_empty & resection_or_ablation & good_spikes;

%% Convert pc and coh to mean across electrodes
pc_mean = cellfun(@(x) nanmean(x,2),pc,'uniformoutput',false);
coh_mean = cellfun(@(x) squeeze(nanmean(x,2)),coh,'uniformoutput',false);
coh_delta_mean = cellfun(@(x) x(:,1),coh_mean,'uniformoutput',false);
coh_theta_mean = cellfun(@(x) x(:,2),coh_mean,'uniformoutput',false);
coh_alpha_mean = cellfun(@(x) x(:,3),coh_mean,'uniformoutput',false);
coh_beta_mean = cellfun(@(x) x(:,4),coh_mean,'uniformoutput',false);
coh_gamma_mean = cellfun(@(x) x(:,5),coh_mean,'uniformoutput',false);

%% Normalize things
rate = cellfun(@(x) (x-nanmean(x))/nanstd(x),rate,'uniformoutput',false);

%% Transform rl
rate_min = 0.1;
rl_txf = cellfun(@(x,y) transform_rl(x,y,rate_min),spike_rl,rate,'uniformoutput',false);
rl = rl_txf;

%% How does rl correlate with spike rate? (a fair amount)
rl_rate_corr = cellfun(@(x,y) corr(x,y,'rows','pairwise'),rl_txf,rate);

%% Make a patient indicator
ps = num2cell((1:npts)');
pid = cellfun(@(x,y) y*ones(length(x),1), rate,ps,'uniformoutput',false);

%% Outcome indicator
good_outcome = complete & outcome_num == 1;
good_outcome_id = num2cell(good_outcome);
gid = cellfun(@(x,y) y*ones(length(x),1), rate,good_outcome_id,'uniformoutput',false);

%% Make table
T = table(cat(1,soz{:}),cat(1,rl{:}),cat(1,rate{:}),cat(1,pid{:}),...
    cat(1,gid{:}),cat(1,pc_mean{:}),...
    cat(1,coh_delta_mean{:}),cat(1,coh_theta_mean{:}),...
    cat(1,coh_alpha_mean{:}),cat(1,coh_beta_mean{:}),cat(1,coh_gamma_mean{:}),...
    'VariableNames',...
    {'SOZ','rl','rate','pid','good_outcome','pc','delta_coh',...
    'theta_coh','alpha_coh','beta_coh','gamma_coh'});


%% Design classifier parameters
params.Nsplits = 1e2;
params.perc_train = 2/3;
params.predictors = {'rate'};
params.response = 'SOZ';
params.pt_id = 'pid';
params.loo = 0;

all_auc = general_logistic_regression(T,params);
nanmean(all_auc)

%{
unpaired_plot(all_auc(complete & outcome_num == 1),all_auc(complete & outcome_num == 0),...
    {'Good outcome','Bad outcome'},'AUC')
%}
end