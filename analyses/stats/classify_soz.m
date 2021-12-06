function classify_soz

%{
mixed effects Logistic regression model
SOZ ~ rate_wake + rate_sleep + rl_wake + rl_sleep + (1|patient)
%}

%% Seed rng (because it's a MC test)
rng(0)
nb = 1e4;

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
rate_sw = out.bin_out.all_elec_rates_sw;
rl_sw = out.bin_out.all_elecs_rl_sw;
soz = out.bin_out.all_is_soz;
locs = out.circ_out.all_locs;
pt_idx = (1:length(rate_sw))';

%% Put data into friendly format
vec_rate
for ip = 1:length(pt_idx)
    
    
end

%% Divide into training and testing set
training = randsample(length(pt_idx),floor(length(pt_idx)/2));
training_idx = ismember(pt_idx,training);
testing_idx = ~ismember(pt_idx,training);





end