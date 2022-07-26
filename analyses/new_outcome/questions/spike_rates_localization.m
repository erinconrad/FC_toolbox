function spike_rates_localization

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
good_spikes = data.good_spikes;
spike_rates = data.all_spikes;
locs = data.all_soz_locs;
lats = data.all_soz_lats;

%% Parse localizations
temporal = contains(locs,'temporal');
oc = strcmp(locs,'other cortex');
multi = contains(locs,'diffuse') | (contains(locs,'multifocal') & ~contains(locs,'temporal'));

%% Take mean spike rate
mean_spike_rates = cellfun(@nanmean,spike_rates);


spikes_temporal = mean_spike_rates(good_spikes & temporal);
spikes_oc = mean_spike_rates(good_spikes & oc);
spikes_multi = mean_spike_rates(good_spikes & multi);

unpaired_plot(spikes_temporal,spikes_oc,{'Temporal','Other cortex'},'Mean spike rate')



end