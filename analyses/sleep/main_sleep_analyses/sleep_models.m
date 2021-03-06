function nmout = sleep_models(just_gray,pre_implant)

%{
This code executes the SOZ localization models for the spikes and sleep
paper. Just_gray indictaes whether to do the analysis on only gray matter
electrodes.
%}

%% Parameters
durations = {1, 5, 10, 30, 60,[]}; % durations to check for duration analysis


%% Locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))

%% Do the stats

% Bootstrap stats
fprintf('\nDoing bootstrap to get model stats\n');
[coeff_names,coeff_stats] = sleep_model_bootstrap_stats(just_gray,pre_implant);

% LOO
fprintf('\nDoing LOO analysis to get individual patient performance\n');
[pt_stats,X,Y,pt_specific,excluded] = sleep_loo(just_gray,pre_implant);

% Get AUC for sleep and wake as a function of duration
fprintf('\nDoing duration analysis\n');
time_aucs = sleep_duration(durations,just_gray,pre_implant);

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