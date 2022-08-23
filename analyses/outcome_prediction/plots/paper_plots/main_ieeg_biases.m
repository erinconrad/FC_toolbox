function main_ieeg_biases(n_test_train_splits,doing_from_github)

% seed random number generator (for consistency in the random testing/training splits)
rng(0) 

% how many testing/training splits
if ~exist('n_test_train_splits','var')
    n_test_train_splits = 1e3;
end

% 1 means doing from github so skip some conceptual plots and Table 1, which require looking at larger datasets
if ~exist('doing_from_github','var')
    doing_from_github = 1;
end

%% Get file locs and set script path
locations = fc_toolbox_locs;
scripts_folder = locations.script_folder;
plot_folder = locations.paper_plot_folder;
addpath(genpath(scripts_folder));

%% Symmetric coverage tests (both atlases)
fprintf('\nDoing symmetric coverage tests...\n');
symmetric_coverage_tests('brainnetome',0)
symmetric_coverage_tests('aal_bernabei',0)

%% Symmetric coverage tests for post-review analysis
% Testing if I get similar result if I subsample electrodes
nb = 1e2;
for ib = 1:nb
   a = symmetric_coverage_tests('brainnetome',1);
   a = rmfield(a,'soz_non_soz_ordered_atlas'); % large, don't need
   nbrain_out(ib) = a;
   a = symmetric_coverage_tests('aal_bernabei',1);
   a = rmfield(a,'soz_non_soz_ordered_atlas'); % large, don't need
   naal_out(ib) = a;
end
save([plot_folder,'nbrain_out.mat'],'nbrain_out')
save([plot_folder,'naal_out.mat'],'naal_out')

%% SOZ classifier
fprintf('\nDoing classifier to predict SOZ vs non-SOZ (takes a while if many training/testing splits)...\n');
simpler_classifier('brainnetome',n_test_train_splits)
simpler_classifier('aal_bernabei',n_test_train_splits)

%% Plots
fprintf('\nGenerating plots...\n');
%do_all_figs(doing_from_github) 

end