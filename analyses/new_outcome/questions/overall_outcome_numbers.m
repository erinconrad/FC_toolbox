function overall_outcome_numbers

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;

 %% get variables of interest
ilae = data.all_two_year_ilae;
engel = data.all_two_year_engel;
surgery = data.all_surgery;
good_spikes = data.good_spikes;

%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);

%% Find those with non-empty outcomes
non_empty_engel = cellfun(@(x) ~isempty(x), engel);
non_empty_ilae = cellfun(@(x) ~isempty(x), ilae);

%% Find those with non empty outcomes, resection or ablation, and good spikes
engel_complete = non_empty_engel & resection_or_ablation & good_spikes;
ilae_complete = non_empty_ilae & resection_or_ablation & good_spikes;

fprintf('\nThere are %d patients with good spikes, resection or ablation, and 2 year ilae scores.\n',sum(ilae_complete));
fprintf('\nThere are %d patients with good spikes, resection or ablation, and 2 year engel scores.\n',sum(engel_complete));

%% Get all unique outcomes in both sets
[engel_cats,~,ic_engel] = unique(engel);
[ilae_cats,~,ic_ilae] = unique(ilae);

%% Make a summary table
engel_counts = cellfun(@(x) sum(strcmp(x,engel)),engel_cats);
ilae_counts = cellfun(@(x) sum(strcmp(x,ilae)),ilae_cats);

end