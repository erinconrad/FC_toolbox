function lr_spike_asymmetry_outcome

which_atlas = 'aal';
which_outcome = 'ilae';
restrict_mt = 0;

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
locs = data.all_locs;
labels = data.all_labels;
anatomy = data.all_anatomy;
soz = data.all_soz_bin;

%% get atlas
switch which_atlas
    case 'aal'
        atlas = data.all_aal;
        atlas_names = data.aal_names;
        mt_names = {'Hippocampus','Amygdala'};
    case 'brainnetome'
        atlas = data.all_brainnetome;
        atlas_names = data.brainnetome_names;
        mt_names = {'Amyg','Hipp'};
    
end

%% Define names corresponding to mesial temporal (not used in paper)
%mt = cellfun(@(x) ismember(x,mt_names),atlas,'uniformoutput',false);

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
complete = non_empty_outcome & resection_or_ablation & good_spikes;

%% Find good and bad outcome
outcome_num = cellfun(@(x) parse_outcome(x,which_outcome),outcome);

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
if restrict_mt
    mt_lr = cellfun(@(x,y) return_bilateral_mesial_temporal_spikes(x,y,which_atlas),...
        spike_rates,atlas,'uniformoutput',false);
    mean_lr_spike_rate = cell2mat(mt_lr);
else
    mean_lr_spike_rate = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'))) nanmean(x(strcmp(y,'R')))],...
    spike_rates,elec_lats,'uniformoutput',false);
    mean_lr_spike_rate = cell2mat(mean_lr_spike_rate);
end

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

% Weighted dispersion of soz
SDSz = cellfun(@(x,y) weighted_standard_distance(x,y),locs,soz);
SDSznorm = SDSz./SDE;

%% Does weighted dispersion correlate with ai?
% yeah pretty negative correlation, implying that a greater asymmetry index
% implies LOWER weighted dispersion (spikes are more focal)
if 0
figure
plot(ai,SD,'o')
end

%% Does spike spatial dispersion correlate with electrode spatial dispersion and SOZ spatial dispersion?
% Yes, so important to control for the others when predicting outcome
if 0 
   figure
   nexttile
   plot(SD,SDE,'o')
    
   nexttile
   plot(SD,SDSz,'o')
end

%% Plot
figure
set(gcf,'position',[440 442 1400 355])
tiledlayout(1,3,'tilespacing','tight','padding','tight')
nexttile
stats = unpaired_plot(ai(outcome_num==1 & complete),...
    ai(outcome_num==0 & complete),{'Good outcome','bad outcome'},'Mean abs diff spike rate');
set(gca,'fontsize',20)
title('Spike asymmetry')

nexttile
stats2 = unpaired_plot(SDnorm(outcome_num==1 & complete),...
    SDnorm(outcome_num==0 & complete),{'Good outcome','bad outcome'},'Spike spatial dispersion');
set(gca,'fontsize',20)
title('Spike spatial spread')

nexttile
stats3 = unpaired_plot(SDSznorm(outcome_num==1 & complete),...
    SDSznorm(outcome_num==0 & complete),{'Good outcome','bad outcome'},'SOZ spatial dispersion');
set(gca,'fontsize',20)
title('SOZ spatial spread')


end