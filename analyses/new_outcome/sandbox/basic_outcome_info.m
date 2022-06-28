function basic_outcome_info

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
data_folder = [results_folder,'analysis/new_outcome/data/'];


% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Load data file
data = load([data_folder,'main_out.mat']);
data = data.out;

%% outcomes
ilae = data.all_two_year_ilae;
engel = data.all_two_year_engel;

%% Get all unique outcomes in both sets
[engel_cats,~,ic_engel] = unique(engel);
[ilae_cats,~,ic_ilae] = unique(ilae);

%% Make a summary table
engel_counts = cellfun(@(x) sum(strcmp(x,engel)),engel_cats);
ilae_counts = cellfun(@(x) sum(strcmp(x,ilae)),ilae_cats);

table(engel_cats,engel_counts)
table(ilae_cats,ilae_counts)

end