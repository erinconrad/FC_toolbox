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
rate = out.bin_out.all_elec_rates;
rl = out.bin_out.all_elecs_rl;
soz_rank_sw_rate = out.bin_out.soz_rank_sw;
soz_rank_sw_rl = out.bin_out.soz_rank_sw_rl;
soz = out.bin_out.all_is_soz;
locs = out.circ_out.all_locs;


end