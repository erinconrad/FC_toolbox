%{
This is the primary script to run all the analyses for the epilepsy
laterality project.
%}

do_full_pipeline = 0; % switch to 1 if you want to run the full analysis (takes several hours)

% End users should not change this
doing_from_github = 1;

%% Get file locs
locations = fc_toolbox_locs;
plot_folder = locations.el_plots_folder;

% Remove the results.html file
if exist([plot_folder,'results.html'],'file') ~= 0
    delete([plot_folder,'results.html'])
end

if exist([plot_folder,'supplemental_results.html'],'file') ~= 0
    delete([plot_folder,'supplemental_results.html'])
end

%% Do model analyses
% This takes the intermediate file containing electrode-level features,
% calculates AI, and then does models to classify SOZ. This takes a while,
% in large part because it needs to run the subsampling analyses (takes
% 6 hours or so). It generates .mat files in the plot folder containing the
% model results. If this is set to 0, then the remaining parts of the
% function use the existing .mat model result files to generate figures. If
% set to 1, it replaces these .mat files with new ones.
if do_full_pipeline
    model_ext_validation_files
end

%% Now do figures

% Figure 1 is methods figure

% Table 1
if ~doing_from_github
    % This requires loading files that have a lot of clinical data and are
    % too large to share
    make_table_1
end

% Table S1 is feature description table

% Figure 2
combined_univariate_fmri_plots

% Figure S1
correlation_figure

% Figure 3 and Fig S2 and S3
model_plots

% Figure 4 and Fig S4 and S5
outcome_plots

close all

