function regional_spikes_outcome

%% Parameters
which_outcome = 'ilae';
which_atlas = 'aal'; 
granularity = 'spatial';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];

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
aal_atlas_names = data.aal_names;
brainnetome_atlas_names = data.brainnetome_names;
aal = data.all_aal;
anatomy = data.all_anatomy;
brainnetome = data.all_brainnetome;
npts = length(good_spikes);
locs = data.all_locs;


switch which_atlas
    case 'aal'
        atlas = aal;
        atlas_names = aal_atlas_names;
    case 'brainnetome'
        atlas = brainnetome;
        atlas_names = brainnetome_atlas_names;
    
end

%% Get outcome
switch which_outcome
    case 'ilae'
        outcome = ilae;
    case 'engel'
        outcome = engel;
end

%% Find good and bad outcome
outcome_num = cellfun(@(x) parse_outcome(x,which_outcome),outcome);
outcome_cat = cell(length(outcome_num),1);
outcome_cat(outcome_num==1) = {'good'};
outcome_cat(outcome_num==0) = {'bad'};

%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);

%% Find those with non-empty outcomes
non_empty = cellfun(@(x) ~isempty(x), outcome);
complete = resection_or_ablation & non_empty & good_spikes;

%% convert spikes to atlas space
spikes_bin = cellfun(@(x,y) bin_univariate_atlas(x,y,atlas_names),...
    spike_rates,atlas,'uniformoutput',false);
spikes_bin = horzcat(spikes_bin{:}); % I still havent removed bad spikes, but they should all be nans

%% Break atlas into categories
broad_regions = localize_regions(atlas_names,which_atlas);
non_empty_broad_regions = (cellfun(@(x) ~isempty(x),broad_regions));
non_empty_names = broad_regions(non_empty_broad_regions);
broad_regions(non_empty_broad_regions) = non_empty_names;
broad_regions(~non_empty_broad_regions) = {''};

%% Further broaden broad regions
final_regions = broad_regions;

switch granularity
    case 'fine'
        % no changes
    case 'lobar'
        % combine MT and NC
        final_regions = cellfun(@(x) strrep(x,'mesial temporal','temporal'),...
            final_regions,'uniformoutput',false);
        final_regions = cellfun(@(x) strrep(x,'temporal neocortical','temporal'),...
            final_regions,'uniformoutput',false);
    case 'hemispheric'
        % combine all left, all right
        right_regions = cellfun(@(x) contains(x,'right'),final_regions);
        left_regions = cellfun(@(x) contains(x,'left'),final_regions);
        final_regions(right_regions) = {'right'};
        final_regions(left_regions) = {'left'};
    case 'individual'
end

if 0
    table(broad_regions,final_regions)
end

%% Get broad regional identities of each electrode
elec_broad_names = cellfun(@(x) elec_broad(x,atlas_names,final_regions),...
    atlas,'uniformoutput',false);

%% Get average spike rates in these broad regions
unique_regions = unique(final_regions(cellfun(@(x) ~isempty(x),final_regions)));
nregions = length(unique_regions);
spikes_broad = nan(nregions,npts);
for i = 1:nregions
    spikes_broad(i,:) = nanmean(spikes_bin(strcmp(final_regions,unique_regions{i}),:),1);
end

%% Get number of elecs per broad region
num_elecs = cellfun(@(x) num_elecs_region(unique_regions,x),...
    elec_broad_names,'uniformoutput',false);
num_elecs = horzcat(num_elecs{:});

%% Define L-R lateralizations
atlas_lat = lateralize_regions_simple(atlas_names);
elec_lats = cellfun(@(x) elec_broad(x,atlas_names,atlas_lat), atlas,'uniformoutput',false);

if 0
    %table(atlas_lat,atlas_names)
    i = 100;
    table(elec_lats{i},anatomy{i})
end

%% Get average left and right spike rates
mean_lr_spike_rate = cellfun(@(x,y) [nanmean(x(strcmp(y,'L'))) nanmean(x(strcmp(y,'R')))],...
    spike_rates,elec_lats,'uniformoutput',false);
mean_lr_spike_rate = cell2mat(mean_lr_spike_rate);



%% Clean the data a little
spikes_broad = spikes_broad';
num_elecs = num_elecs';

complete = complete & sum(num_elecs,2) ~= 0;

% remove bad patients
spikes_broad(~complete,:) = [];
num_elecs(~complete,:) = [];
outcome_num(~complete) = [];
outcome_cat(~complete) = [];
spike_rates(~complete) = [];
mean_lr_spike_rate(~complete,:) = [];
locs(~complete,:) = [];

% nan spikes are exactly those without elecs
assert(isequal(find(isnan(spikes_broad)),find(num_elecs==0)))

% make nans zero
%spikes_broad(isnan(spikes_broad)) = 0;


%% Define some regional measure of asymmetry
% How localized are the spikes to one of the N regions
if strcmp(granularity,'spatial')
    %% Also get weighted dispersion of spikes
    SD = cellfun(@(x,y) weighted_standard_distance(x,y),locs,spike_rates);
    SD(SD>1e10) = nan;

    % Also weighted dispersion of electrodes
    SDE = cellfun(@(x) weighted_standard_distance(x,[]),locs);
    SDE(SDE>1e10) = nan;
    SDnorm = SD./SDE;
    sp_asymm = SD;
    elec_asymmetry = SDE;
    AI = SDnorm;
elseif strcmp(granularity,'individual')
    AI = cellfun(@regional_spike_asymmetry,spike_rates);
elseif strcmp(granularity,'broad')
    spikes_broad_cell = (num2cell(spikes_broad',1))';
    AI = cellfun(@regional_spike_asymmetry,spikes_broad_cell);
elseif strcmp(granularity,'hemispheric')
    %{
    hemispheric_broad(:,1) = nanmean(spikes_broad(:,contains(unique_regions,'left')),2);
    hemispheric_broad(:,2) = nanmean(spikes_broad(:,contains(unique_regions,'right')),2);
    hemispheric_broad_cell = (num2cell(hemispheric_broad',1))';
    %}
    hemispheric_broad_cell = (num2cell(mean_lr_spike_rate',1))';
    AI = cellfun(@regional_spike_asymmetry,hemispheric_broad_cell);
    
    % Double check that my regional asymmetry gives same answer as AI for
    % hemispheric condition
    %{
    L = hemispheric_broad(:,1);
    R = hemispheric_broad(:,2);
    alt_ai = abs(L-R)./(L+R);
    assert(nansum(abs(AI-alt_ai))<1e-5)
    %}
end

if 0
%% Compare AI for good and bad outcome
unpaired_plot(AI(outcome_num==1),AI(outcome_num==0),{'Good outcome','Bad outcome'},'AI')
end

%% Now define broad spike rates to be zero for PCA
spikes_broad(isnan(spikes_broad)) = 0;

%% Normalize spikes broad
spikes_broad = spikes_broad./sum(spikes_broad,2);

%% Construct table
full = [outcome_num,spikes_broad,num_elecs];
predictorNames = [(cellfun(@(x) sprintf('%s_spikes',x),unique_regions,'uniformoutput',false))',...
    (cellfun(@(x) sprintf('%s_elecs',x),unique_regions,'uniformoutput',false))'];
responseName = 'Outcome';
T = array2table(full,'variablenames',[responseName,predictorNames]);

null = [outcome_num,num_elecs];
predictorNamesNull = [(cellfun(@(x) sprintf('%s_elecs',x),unique_regions,'uniformoutput',false))'];
responseName = 'Outcome';
TNull = array2table(null,'variablenames',[responseName,predictorNamesNull]);

%% Train/test classifier
[trueClass,predClass,AUC] = more_general_train_test(T,1e2,2/3,@log_regression_PCA,predictorNames,responseName);
[trueClassNull,predClassNull,AUCNull] = more_general_train_test(TNull,1e2,2/3,@log_regression_PCA,predictorNamesNull,responseName);

mean(AUC)
mean(AUCNull)
end