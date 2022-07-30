function spike_correspondence_atl_outcome

which_atlas = 'aal';
which_outcome = 'ilae';
atl_or_mt_ablation = 'atl';

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
npts = length(good_spikes);
spike_rates = data.all_spikes;
resec_lat = data.all_resec_lat;
resec_loc = data.all_resec_loc;
ablate_lat = data.all_ablate_lat;
ablate_loc = data.all_ablate_loc;
aal_atlas_names = data.aal_names;
brainnetome_atlas_names = data.brainnetome_names;
aal = data.all_aal;
brainnetome = data.all_brainnetome;

switch atl_or_mt_ablation
    case 'atl'
        which_localization = 'broad';
        surg_loc = resec_loc;
        surg_lat = resec_lat;
        loc_of_interest = 'ATL';
        lobe = 'temporal';
    case 'mt_ablation'
        which_localization = 'fine';
        surg_loc = ablate_loc;
        surg_lat = ablate_lat;
        loc_of_interest = 'mesial temporal';
        lobe = 'mesial temporal';
end

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


%% Parse surgery
resection_or_ablation = cellfun(@(x) ...
    contains(x,'resection','ignorecase',true) | contains(x,'ablation','ignorecase',true),...
    surgery);

%% Find those with non-empty outcomes
non_empty_outcome = cellfun(@(x) ~isempty(x), outcome);

%% Define complete
complete = non_empty_outcome & resection_or_ablation;

%% convert spikes to atlas space
spikes_bin = cellfun(@(x,y) bin_univariate_atlas(x,y,atlas_names),...
    spike_rates,atlas,'uniformoutput',false);
spikes_bin = horzcat(spikes_bin{:}); % I still havent removed bad spikes, but they should all be nans

%% Break atlas into categories
broad_regions = localize_regions(atlas_names,which_atlas);
non_empty_broad_regions = (cellfun(@(x) ~isempty(x),broad_regions));
non_empty_names = broad_regions(non_empty_broad_regions);
if strcmp(which_localization,'broad')
    non_empty_names = cellfun(@(x) strrep(x,'mesial temporal','temporal'),...
        non_empty_names,'uniformoutput',false);
    non_empty_names = cellfun(@(x) strrep(x,'temporal neocortical','temporal'),...
        non_empty_names,'uniformoutput',false);
end
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

% Get percentage of elecs in that region
prop_elecs = num_elecs./nansum(num_elecs,1)*100;

%% Also get percentage of spikes in that region
perc_spikes_broad = spikes_broad./nansum(spikes_broad,1)*100;

%% Define proportion of elecs in temporal lobe ipsilateral to ATL
% And proportion of spikes in temporal lobe ipsilateral to ATL
prop_elecs_atl = nan(npts,1);
prop_spikes_atl = nan(npts,1);
for ip = 1:npts
    if ~strcmp(surg_loc{ip},loc_of_interest)
        continue
    end
    
    if strcmp(surg_lat{ip},'right')
        reg = sprintf('right %s',lobe);
    elseif strcmp(surg_lat{ip},'left')
        reg = sprintf('left %s',lobe);
    else
        continue
    end
    
    prop_elecs_atl(ip) = prop_elecs(strcmp(unique_regions,reg),ip);
    prop_spikes_atl(ip) = perc_spikes_broad(strcmp(unique_regions,reg),ip);
end


%% Plot
stats = unpaired_plot(prop_spikes_atl(outcome_num==1 & complete),...
    prop_spikes_atl(outcome_num==0 & complete),{'Good outcome','bad outcome'},'Proportion of spikes in ipsilateral region');
set(gca,'fontsize',20)

end