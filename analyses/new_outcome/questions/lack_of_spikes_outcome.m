function lack_of_spikes_outcome

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
spike_rates = data.all_spikes;

%% Are all spikes nans if not good? - Yes
assert(isequal(cellfun(@(x) all(isnan(x)),spike_rates),~good_spikes))

%% Take mean spike rate
mean_spike_rates = cellfun(@nanmean,spike_rates);

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

%% Find good and bad ilae outcome
ilae_num = cellfun(@(x) parse_outcome(x,'ilae'),ilae);

%% Compare spike rates in good versus bad outcome
spikes_good = mean_spike_rates(ilae_complete & ilae_num == 1);
spikes_bad = mean_spike_rates(ilae_complete & ilae_num == 0);

%% Plot
unpaired_plot(spikes_good,spikes_bad,{'Good outcome','bad outcome'},'Mean spike rate')

end