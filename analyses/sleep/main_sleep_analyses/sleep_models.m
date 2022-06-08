function nmout = sleep_models(just_gray)

%% Parameters
durations = {1, 5, 10, 30, 60,[]};


%% Locations
locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/epilepsia/'];
out_folder1 = [script_folder,'analyses/sleep/data/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

%% Load out file and get roc stuff
out = load([out_folder1,'out.mat']);
out = out.out;



%% Do the stats

% Bootstrap stats
fprintf('\nDoing bootstrap to get model stats\n');
[coeff_names,coeff_stats] = sleep_model_bootstrap_stats(just_gray);

% LOO
fprintf('\nDoing LOO analysis to get individual patient performance\n');
[pt_stats,X,Y,pt_specific,excluded] = sleep_loo(just_gray);

% Get AUC for sleep and wake as a function of duration
fprintf('\nDoing duration analysis\n');
time_aucs = sleep_duration(durations,just_gray);

nmout.coeff_stats = coeff_stats;
nmout.coeff_names = coeff_names;
nmout.pt_stats = pt_stats;
nmout.X = X;
nmout.Y = Y;
nmout.pt_specific = pt_specific;
nmout.time_aucs = time_aucs;
nmout.excluded = excluded;
nmout.durations = durations;



        
    


end