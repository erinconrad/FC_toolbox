function do_all_results

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
plot_folder = [results_folder,'analysis/new_outcome/plots/'];

% Remove the results.html file
if exist([plot_folder,'results.html'],'file') ~= 0
    delete([plot_folder,'results.html'])
end

if exist([plot_folder,'supplemental_results.html'],'file') ~= 0
    delete([plot_folder,'supplemental_results.html'])
end

% Figure 1 is methods figure

% Table 1
make_table_1

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

end