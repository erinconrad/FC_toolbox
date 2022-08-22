function spikes_predict_bilaterality

which_atlas = 'brainnetome';

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
locs = data.all_locs;
labels = data.all_labels;

anatomy = data.all_anatomy;
soz = data.all_soz_bin;
soz_lats = data.all_soz_lats;

%% get atlas
switch which_atlas
    case 'aal'
        atlas = data.all_aal;
        atlas_names = data.aal_names;
    case 'brainnetome'
        atlas = data.all_brainnetome;
        atlas_names = data.brainnetome_names;
    
end

%% Define bilateral patients
bilat = strcmp(soz_lats,'bilateral') | strcmp(soz_lats,'diffuse');
unilat = strcmp(soz_lats,'left') | strcmp(soz_lats,'right');


%% Are all spikes nans if not good? - Yes
assert(isequal(cellfun(@(x) all(isnan(x)),spike_rates),~good_spikes))

%% Define L-R lateralizations
if strcmp(which_atlas,'labels')
    elec_lats = cellfun(@(x,y) label_and_anatomy_lat_determination(x,y),labels,anatomy,'uniformoutput',false);
else
    atlas_lat = lateralize_regions_simple(atlas_names);
    elec_lats = cellfun(@(x) elec_broad(x,atlas_names,atlas_lat), atlas,'uniformoutput',false);
end

%% Get average left and right spike rates
mean_lr_spike_rate = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'))) nanmean(x(strcmp(y,'R')))],...
    spike_rates,elec_lats,'uniformoutput',false);
mean_lr_spike_rate = cell2mat(mean_lr_spike_rate);


%% Get two different asymmetric indices
L = mean_lr_spike_rate(:,1);
R = mean_lr_spike_rate(:,2);
ai = max([L./(L+R),R./(L+R)],[],2);
alt_ai = abs(L-R)./(L+R);

%% Also get weighted dispersion of spikes
SD = cellfun(@(x,y) weighted_standard_distance(x,y),locs,spike_rates);
SD(SD>1e10) = nan;

% Also weighted dispersion of electrodes
SDE = cellfun(@(x) weighted_standard_distance(x,[]),locs);
SDE(SDE>1e10) = nan;
SDnorm = SD./SDE;

unpaired_plot(alt_ai(unilat),alt_ai(bilat),{'Unilateral','Bilateral'},'Asymmetry index')


end