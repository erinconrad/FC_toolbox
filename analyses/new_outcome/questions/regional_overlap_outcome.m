function regional_overlap_outcome

%{
This attempts to answer whether a higher proportion of spikes coming from
the region of the soz predicts a better outcome.

1. Get SOZ region
2. Calc prop of spikes in that region/prop of elecs in that region
3. See if this is higher in patients with good outcome

What to do with broader SOZs?
- If totally diffuse, then throw the patient out
- If left diffuse or right diffuse, could include and then just compare
left to right spikes

ALSO, think about how granular I want to be. Could just be lobe. Could
even just start with hemisphere...

%}

%% Parameters
which_outcome = 'ilae';
which_atlas = 'aal'; 
granularity = 'fine';

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
soz = data.all_soz_bin;
aal_atlas_names = data.aal_names;
brainnetome_atlas_names = data.brainnetome_names;
aal = data.all_aal;
brainnetome = data.all_brainnetome;
npts = length(good_spikes);
soz_locs = data.all_soz_locs;
soz_lats = data.all_soz_lats;


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

%% convert spikes to atlas space
spikes_bin = cellfun(@(x,y) bin_univariate_atlas(x,y,atlas_names),...
    spike_rates,atlas,'uniformoutput',false);
spikes_bin = horzcat(spikes_bin{:}); % I still havent removed bad spikes, but they should all be nans

%% Break atlas into categories
broad_regions = localize_regions(atlas_names,which_atlas);
non_empty_broad_regions = (cellfun(@(x) ~isempty(x),broad_regions));
non_empty_names = broad_regions(non_empty_broad_regions);
broad_regions(non_empty_broad_regions) = non_empty_names;


%% Get broad regional identities of each electrode
elec_broad_names = cellfun(@(x) elec_broad(x,atlas_names,broad_regions),...
    atlas,'uniformoutput',false);

%% Get average spike rates in these broad regions
unique_regions = unique(broad_regions(cellfun(@(x) ~isempty(x),broad_regions)));
nregions = length(unique_regions);
spikes_broad = nan(nregions,npts);
for i = 1:nregions
    spikes_broad(i,:) = nanmean(spikes_bin(strcmp(broad_regions,unique_regions{i}),:),1);
end


%% Get number of elecs per broad region
num_elecs = cellfun(@(x) num_elecs_region(unique_regions,x),...
    elec_broad_names,'uniformoutput',false);
num_elecs = horzcat(num_elecs{:});


%% Homogenize and combine soz locs and lats
comb = cellfun(@(x,y) homogenize_soz_locs_lats(x,y,'fine'),soz_locs,soz_lats,'uniformoutput',false);


%% Calculate the proportion of spikes and electrodes in the SOZ region
prop_spikes_soz = nan(npts,1);
prop_elecs_soz = nan(npts,1);
for ip = 1:npts
    % Get soz
    curr_soz = comb{ip};
    
    % Get number of spikes and elecs in each region for this patient
    curr_spikes = spikes_broad(:,ip);
    curr_elecs = num_elecs(:,ip);
    
    % make nans zeros
    curr_spikes(isnan(curr_spikes)) = 0;
    curr_elecs(isnan(curr_elecs)) = 0;
    
    %% define SOZ regions
    % Here I can decide how broad to be
    switch granularity
        case 'fine'
            if strcmp(curr_soz,' ') || strcmp(curr_soz,'bilateral/diffuse')
                continue;
            elseif strcmp(curr_soz,'right multifocal/diffuse')
                soz_regions = {'right mesial temporal','right other cortex','right temporal neocortical'};
            elseif strcmp(curr_soz,'left multifocal/diffuse')
                soz_regions = {'left mesial temporal','left other cortex','left temporal neocortical'};
            else
                soz_regions = curr_soz;
            end
            
        case 'lobar'
            if strcmp(curr_soz,' ') || strcmp(curr_soz,'bilateral/diffuse')
                continue;
            elseif strcmp(curr_soz,'right multifocal/diffuse')
                soz_regions = {'right mesial temporal','right other cortex','right temporal neocortical'};
            elseif strcmp(curr_soz,'left multifocal/diffuse')
                soz_regions = {'left mesial temporal','left other cortex','left temporal neocortical'};
            elseif contains(curr_soz,'left') && contains(curr_soz,'temporal')
                soz_regions = {'left mesial temporal','left temporal neocortical'};
            elseif contains(curr_soz,'right') && contains(curr_soz,'temporal')
                soz_regions = {'right mesial temporal','right temporal neocortical'};
            else
                assert(strcmp(curr_soz,'left other cortex') || strcmp(curr_soz,'right other cortex'))
                soz_regions = curr_soz;
            end
            
        case 'hemispheric'
            if strcmp(curr_soz,' ') || strcmp(curr_soz,'bilateral/diffuse')
                continue;
            elseif contains(curr_soz,'left')
                soz_regions = {'left mesial temporal','left other cortex','left temporal neocortical'};
            else
                assert(contains(curr_soz,'right'))
                soz_regions = {'right mesial temporal','right other cortex','right temporal neocortical'};
            end
            
    end
    
    % find matches
    matching_regions = ismember(unique_regions,soz_regions);
    
    % Get proportion of spikes and elecs that are in the matching regions
    curr_prop_spikes = sum(curr_spikes(matching_regions))./sum(curr_spikes);
    curr_prop_elecs = sum(curr_elecs(matching_regions))./sum(curr_elecs);
    
    prop_spikes_soz(ip) = curr_prop_spikes;
    prop_elecs_soz(ip) = curr_prop_elecs;
    
   
end

% Get ratio of proportion of spikes to proportion of elecs
rat_prop = prop_spikes_soz./prop_elecs_soz;

%% get the patients I want
include = resection_or_ablation & non_empty & good_spikes;

%% Plot comparison of ratio of prop spikes to prop elecs for good vs bad outcome
unpaired_plot(rat_prop(include&outcome_num==1),rat_prop(include&outcome_num==0),...
    {'Good outcome','Bad outcome'},'Ratio of prop spikes in soz to prop elecs in soz')
set(gca,'fontsize',15)

end