function main_ieeg_biases(n_test_train_splits,doing_from_github)

% seed random number generator (for consistency in the random testing/training splits)
rng(0) 

% how many testing/training splits
if ~exist('n_test_train_splits','var')
    n_test_train_splits = 1e3;
end

% 1 means doing from github so skip some conceptual plots and Table 1 that require looking at larger datasets
if ~exist('doing_from_github','var')
    doing_from_github = 1;
end

%% Get file locs and set script path
locations = fc_toolbox_locs;
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Symmetric coverage tests (both atlases)
fprintf('\nDoing symmetric coverage tests...\n');
symmetric_coverage_tests('brainnetome')
symmetric_coverage_tests('aal_bernabei')

%% Spike-FC correlation
fprintf('\nDoing spike analyses...\n');
spike_fc_correlation

%% Density model
fprintf('\nFinding optimal search radius for density model...\n');
erin_dens_model

%% SOZ classifier
fprintf('\nDoing classifier to predict SOZ vs non-SOZ (takes a while if many training/testing splits)...\n');
soz_classifier('brainnetome',n_test_train_splits)
soz_classifier('aal_bernabei',n_test_train_splits)

%% Plots
fprintf('\nGenerating plots...\n');
do_all_figs(doing_from_github) 

end