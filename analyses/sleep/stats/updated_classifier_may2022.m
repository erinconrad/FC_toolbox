function updated_classifier_may2022(leave_out)

%% Parameters
prop_train = 2/3;
do_norm = 1;
randomize_soz = 0; % negative control. If I shuffle the SOZ, I should not exceed chance AUC


locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
script_folder = locations.script_folder;
results_folder = [locations.main_folder,'results/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];


%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;

%% Get stuff
rate_sw = out.bin_out.all_elecs_rates_sw; %rates wake and sleep
soz = out.bin_out.all_is_soz;
pt_idx = (1:length(rate_sw))';
rate_post = out.sz_out.post_ictal_rates;

%% Put data into friendly format
vec_rate_sleep = [];
vec_rate_wake = [];
vec_rate_post = [];
vec_pt_idx = [];
vec_soz = [];
for ip = 1:length(pt_idx)
    
    curr_rate_sw = rate_sw{ip};
    curr_rate_post = rate_post{ip};
    
    if do_norm
        % normalize the rate across electrodes (this is so that patients
        % with higher spike rates in general will get the same weight as
        % patients with lower spike rates)
        curr_rate_sw = (curr_rate_sw - nanmean(curr_rate_sw,1))./...
            nanstd(curr_rate_sw,[],1);        
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


end