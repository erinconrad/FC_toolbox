function lr_spike_asymmetry_outcome

which_atlas = 'aal';
which_outcome = 'ilae';

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

%% get atlas
switch which_atlas
    case 'aal'
        atlas = data.all_aal;
        atlas_names = data.aal_names;
    case 'brainnetome'
        atlas = data.all_brainnetome;
        atlas_names = data.brainnetome_names;
    
end

%% Get outcome
switch which_outcome
    case 'ilae'
        outcome = ilae;
    case 'engel'
        outcome = engel;
end

%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);

%% Find those with non-empty outcomes
non_empty_outcome = cellfun(@(x) ~isempty(x), outcome);

%% Define complete
complete = non_empty_outcome & resection_or_ablation;

%% Find good and bad outcome
outcome_num = cellfun(@(x) parse_outcome(x,which_outcome),outcome);

%% Are all spikes nans if not good? - Yes
assert(isequal(cellfun(@(x) all(isnan(x)),spike_rates),~good_spikes))

%% Define L-R lateralizations
atlas_lat = lateralize_regions_simple(atlas_names);
elec_lats = cellfun(@(x) elec_broad(x,atlas_names,atlas_lat), atlas,'uniformoutput',false);

%% Get average left and right spike rates
mean_lr_spike_rate = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'))) nanmean(x(strcmp(y,'R')))],...
    spike_rates,elec_lats,'uniformoutput',false);
mean_lr_spike_rate = cell2mat(mean_lr_spike_rate);


%% Get two different asymmetric indices
L = mean_lr_spike_rate(:,1);
R = mean_lr_spike_rate(:,2);
ai = max([L./(L+R),R./(L+R)],[],2);
alt_ai = abs(L-R)./(L+R);


%% Plot
stats = unpaired_plot(ai(outcome_num==1 & complete),...
    ai(outcome_num==0 & complete),{'Good outcome','bad outcome'},'Mean abs diff spike rate');
set(gca,'fontsize',20)


end