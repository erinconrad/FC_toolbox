function spike_loc_correspondence

% Seed RNG
rng(0)

%% Parameters
which_atlas = 'aal'; 
which_localization = 'broad';
N = 1e2;
split = 2/3;

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
aal_atlas_names = data.aal_names;
brainnetome_atlas_names = data.brainnetome_names;
aal = data.all_aal;
brainnetome = data.all_brainnetome;
npts = length(good_spikes);
soz_locs = data.all_soz_locs;
soz_lats = data.all_soz_lats;
anatomy = data.all_anatomy;

switch which_atlas
    case 'aal'
        atlas = aal;
        atlas_names = aal_atlas_names;
    case 'brainnetome'
        atlas = brainnetome;
        atlas_names = brainnetome_atlas_names;
    
end

%% convert spikes to atlas space
spikes_bin = cellfun(@(x,y) bin_univariate_atlas(x,y,atlas_names),...
    spike_rates,atlas,'uniformoutput',false);
spikes_bin = horzcat(spikes_bin{:}); % I still havent removed bad spikes, but they should all be nans


%% Confirm patients with bad spikes have all nans for spike rates
assert(isequal(~good_spikes,cellfun(@(x) all(isnan(x)),spike_rates)))

%% Look at the atlas
if 0
figure
turn_nans_gray(spikes_bin)
yticks(1:size(spikes_bin,1))
yticklabels(atlas_names)
end

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

%% Compare atlas names to regions
if 0
    table(broad_regions,atlas_names)
    table(elec_broad_names{1},anatomy{1})
end

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

%% SHow it
if 0
    figure
    turn_nans_gray(num_elecs)
    yticks(1:size(num_elecs,1))
    yticklabels(unique_regions)
end

%% Also get percentage of spikes in that region
perc_spikes_broad = spikes_broad./nansum(spikes_broad,1)*100;

%% Look at it
if 0
figure
turn_nans_gray(spikes_broad)
yticks(1:size(spikes_broad,1))
yticklabels(unique_regions)    
end

%% Homogenize and combine soz locs and lats
comb = cellfun(@(x,y) homogenize_soz_locs_lats(x,y,which_localization),soz_locs,soz_lats,'uniformoutput',false);

%% Build matrix of spike rates by location and SOZ localization
soz_names = unique(comb);
soz_names(strcmp(soz_names,' ')) = [];
nsozs = length(soz_names);

spikes_abs = nan(nsozs,nregions);
spikes_prop = nan(nsozs,nregions);
elecs_prop = nan(nsozs,nregions);
for is = 1:nsozs
    for ir = 1:nregions
        spikes_abs(is,ir) = nanmean(spikes_broad(ir,strcmp(soz_names{is},comb)));
        spikes_prop(is,ir) = nanmean(perc_spikes_broad(ir,strcmp(soz_names{is},comb)));
        elecs_prop(is,ir) = nanmean(prop_elecs(ir,strcmp(soz_names{is},comb)));
    end
end

%% show it
if 0
figure
set(gcf,'position',[10 10 1400 400])
tiledlayout(1,2,'tilespacing','tight')

nexttile
turn_nans_gray(spikes_prop)
yticks(1:nsozs)
yticklabels(soz_names)
xticks(1:nregions)
xticklabels(unique_regions)
title('% Spikes in region')
set(gca,'fontsize',15)
c = colorbar;
ylabel(c,'%','fontsize',15)
ylabel('SOZ localization')
xlabel('Spike region')

nexttile
turn_nans_gray(elecs_prop)
yticks(1:nsozs)
yticklabels(soz_names)
xticks(1:nregions)
xticklabels(unique_regions)
title('% Electrode contacts in region')
set(gca,'fontsize',15)
c = colorbar;
ylabel(c,'%','fontsize',15)
ylabel('SOZ localization')
xlabel('Spike region')
    
end

%% Remove patients with all missing data

% Find patients without atlas localizations
no_atlas = cellfun(@(x) all(cellfun(@(y) isempty(y),x)),atlas);

% Find patients with no soz
no_soz = strcmp(comb,' ');

% Am I correct that those with all nans for spikes are those with no atlas
% or those with bad spikes? - YES!
assert(isequal(no_atlas | ~good_spikes,all(isnan(perc_spikes_broad),1)'))

% prep those to remove
to_remove = no_atlas | ~good_spikes | no_soz;

% Remove the patients
perc_spikes_broad(:,to_remove) = [];
comb(to_remove) = [];
spikes_broad(:,to_remove) = [];
prop_elecs(:,to_remove) = [];

% confirm no more missing data
assert(sum(all(isnan(perc_spikes_broad),1)') == 0)

%% Make spikes_broad that are nans zeros
% This introduces a bias wherein
perc_spikes_broad(isnan(perc_spikes_broad)) = 0;


%% Put data into table
% Build variable names
elec_num_vars = cellfun(@(x) sprintf('%s proportion elecs',x),unique_regions,'uniformoutput',false);
spike_vars = cellfun(@(x) sprintf('%s proportion spikes',x),unique_regions,'uniformoutput',false);
comb_var = 'SOZ';

T = table(comb,prop_elecs',perc_spikes_broad','variablenames',{comb_var,'Var1','Var2'});
T = splitvars(T,{'Var1','Var2'},'newVariableNames',{elec_num_vars',spike_vars'});

%% Train and test models
% Null model: only takes into account number of electrodes in each location
[CNull,ANull] = train_test(T,N,split,@sozTree,'null');
% Full model: only takes into account number of electrodes in each location
[CFull,AFull] = train_test(T,N,split,@sozTree,'full');

%% Show confusion matrices
figure
set(gcf,'position',[10 10 1400 400])
tiledlayout(1,2,'tilespacing','tight','padding','tight')

% Null
nexttile
meanCNull = mean(CNull,3);
show_confusion(meanCNull,soz_names,...
    'Predicted SOZ','True SOZ','Null model',ANull);

% Full
nexttile
meanCFull = mean(CFull,3);
show_confusion(meanCFull,soz_names,...
    'Predicted SOZ','True SOZ','Full model',AFull);

%% Bootstrap CI?
out = bootstrap_ci_and_p(ANull,AFull);

end