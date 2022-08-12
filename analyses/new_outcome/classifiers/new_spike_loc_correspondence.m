function new_spike_loc_correspondence

%% Parameters
do_pca = 0;
rm_bad_outcome_focal = 0;
which_outcome = 'ilae';
which_atlas = 'aal'; 
predictor_granularity = 'broad';
response_granularity = 'broad';
N = 1e2;
split = 2/3;

%% Seed RNG
rng(0)


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
aal_atlas_names = data.aal_names;
brainnetome_atlas_names = data.brainnetome_names;
aal = data.all_aal;
brainnetome = data.all_brainnetome;
npts = length(good_spikes);
soz_locs = data.all_soz_locs;
soz_lats = data.all_soz_lats;
anatomy = data.all_anatomy;
ilae = data.all_two_year_ilae;
engel = data.all_two_year_engel;
surgery = data.all_surgery;

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
nregions = length(atlas_names);
spikes_bin = cellfun(@(x,y) bin_univariate_atlas(x,y,atlas_names),...
    spike_rates,atlas,'uniformoutput',false);
spikes_bin = horzcat(spikes_bin{:}); % I still havent removed bad spikes, but they should all be nans

%% Same with numbers of elecs
nelecs_bin = nan(nregions,npts);
for ip = 1:npts
    curr_elecs_atlas = atlas{ip};
    for ir = 1:length(atlas_names)
        nelecs_bin(ir,ip) = sum(strcmp(curr_elecs_atlas,atlas_names{ir}));
    end
end


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
if strcmp(predictor_granularity,'broad')
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

%% Handy visualization technique: show all electrodes and their broad regional localization, along with clinical anatomy
if 0
    cmap = colormap(lines(nregions));
    cmap = [0.7 0.7 0.7;cmap];
    i = 115;
    curr_elecs = elec_broad_names{i};
    for j = 1:length(curr_elecs)
        if isempty(curr_elecs{j})
            curr_elecs{j} = '';
        end
    end
    [ia,ib] = ismember(curr_elecs,unique_regions);
    figure
    set(gcf,'position',[10 10 1300 900])
    scatter3(locs{i}(:,1),locs{i}(:,2),locs{i}(:,3),150,ib,'filled');
    hold on
    text(locs{i}(:,1),locs{i}(:,2),locs{i}(:,3),anatomy{i});
    colormap(cmap)
    c = colorbar('ticks',0:nregions,'ticklabels',['none';unique_regions]);
    
end

%% Homogenize and combine soz locs and lats
comb = cellfun(@(x,y) homogenize_soz_locs_lats(x,y,response_granularity),soz_locs,soz_lats,'uniformoutput',false);

%% Get counts in each region


%% Build matrix of spike rates by location and SOZ localization
soz_names = unique(comb);
soz_names(strcmp(soz_names,' ')) = [];
nsozs = length(soz_names);

spikes_abs = nan(nsozs,nregions);
spikes_prop = nan(nsozs,nregions);
elecs_prop = nan(nsozs,nregions);
num_soz = nan(nsozs,1);
for is = 1:nsozs
    num_soz(is) = sum(strcmp(soz_names{is},comb));
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

%% Mark patients with "focal" SOZ designations who have poor outcome or no surgery
non_diffuse = cellfun(@(x) ~contains(x,'diffuse'),comb);
diffuse = cellfun(@(x) contains(x,'diffuse'),comb);

focal_no_surg_or_poor_outcome = non_diffuse & (~resection_or_ablation | ~non_empty | outcome_num == 0);

%% Remove patients with all missing data

% Find patients without atlas localizations
no_atlas = cellfun(@(x) all(cellfun(@(y) isempty(y),x)),atlas);

% Find patients with no soz
no_soz = strcmp(comb,' ');

% Am I correct that those with all nans for spikes are those with no atlas
% or those with bad spikes? - YES!
assert(isequal(no_atlas | ~good_spikes,all(isnan(perc_spikes_broad),1)'))

% prep those to remove
if rm_bad_outcome_focal
    to_remove = no_atlas | ~good_spikes | no_soz | focal_no_surg_or_poor_outcome;
else
    to_remove = no_atlas | ~good_spikes | no_soz;
end


%% Make spikes_broad that are nans zeros
% This introduces a bias wherein there will be no spikes where clinicians
% don't place electrodes. I am controlling for this bias by including
% electrode locations as a null model.
perc_spikes_broad(isnan(perc_spikes_broad)) = 0;
spikes_broad(isnan(spikes_broad)) = 0;

%% Normalize for pca
num_elecs_norm = (num_elecs - mean(num_elecs,1))./nanstd(num_elecs,[],1);
spikes_broad_norm = (spikes_broad - mean(spikes_broad,1))./nanstd(spikes_broad,[],1);

%% Put data into table
% Build variable names
elec_num_vars = cellfun(@(x) sprintf('%s proportion elecs',x),unique_regions,'uniformoutput',false);
spike_vars = cellfun(@(x) sprintf('%s proportion spikes',x),unique_regions,'uniformoutput',false);
comb_var = 'SOZ';

T = table(comb,prop_elecs',perc_spikes_broad','variablenames',{comb_var,'Var1','Var2'});
%T = table(comb,spikes_broad_norm',num_elecs_norm','variablenames',{comb_var,'Var1','Var2'});
T = splitvars(T,{'Var1','Var2'},'newVariableNames',{elec_num_vars',spike_vars'});



%% Also prepare table for original atlas
spikes_bin(isnan(spikes_bin)) = 0;
nelecs_bin(isnan(nelecs_bin)) = 0;
spikes_bin = spikes_bin./sum(spikes_bin,1)*100;
nelecs_bin = spikes_bin./sum(nelecs_bin,1)*100;

elec_num_vars = cellfun(@(x) sprintf('%s proportion elecs',x),atlas_names,'uniformoutput',false);
spike_vars = cellfun(@(x) sprintf('%s proportion spikes',x),atlas_names,'uniformoutput',false);
To = table(comb,spikes_bin',nelecs_bin','variablenames',{comb_var,'Var1','Var2'});
To = splitvars(To,{'Var1','Var2'},'newVariableNames',{elec_num_vars',spike_vars'});

if strcmp(predictor_granularity,'atlas')
    T = To;
    
end

%% Remove things to remove
T_rm = T; T_rm(to_remove,:) = [];

%% Train and test models
if 1
    % Null model: only takes into account number of electrodes in each location
    fprintf('\nTraining null model\n');
    if do_pca
        [CNull,ANull,catsN] = train_test(T_rm,N,split,@sozTreePCA,'null');
    else
        [CNull,ANull,catsN] = train_test(T_rm,N,split,@sozTree,'null');
    end
    fprintf('\nTraining full model\n');
    % Full model: only takes into account number of electrodes in each location
    if do_pca
        [CFull,AFull,catsF] = train_test(T_rm,N,split,@sozTreePCA,'full');
    else
        [CFull,AFull,catsF] = train_test(T_rm,N,split,@sozTree,'full');
    end

    %% Show confusion matrices
    assert(isequal(catsN,catsF))
    assert(isequal(soz_names,catsN))

    figure
    set(gcf,'position',[10 10 1400 400])
    tiledlayout(1,2,'tilespacing','tight','padding','tight')

    % Null
    nexttile
    show_confusion(CNull,soz_names,...
        'Predicted SOZ','True SOZ','Null model',ANull);

    % Full
    nexttile
    show_confusion(CFull,soz_names,...
        'Predicted SOZ','True SOZ','Full model',AFull);
    


    %% Bootstrap CI?
    out = bootstrap_ci_and_p(ANull,AFull);
end

%% Next question - does the regional localization predict outcome for those who had surgery
surg_non_empty = resection_or_ablation & non_empty;
T_out = addvars(T,outcome_num,'NewVariableNames','outcome');
to_remove = ~surg_non_empty | no_atlas | ~good_spikes;
T_out_rm = T_out; T_out_rm(to_remove,:) = [];

[OCNull,OANull] = train_test(T_out_rm,N,split,@outcome_logistic_regression,'null');
fprintf('\nTraining full model\n');
% Full model: only takes into account number of electrodes in each location
[OCFull,OAFull] = train_test(T_out_rm,N,split,@outcome_logistic_regression,'full');

%% Show confusion matrices

figure
set(gcf,'position',[10 10 1400 400])
tiledlayout(1,2,'tilespacing','tight','padding','tight')

% Null
nexttile
show_confusion(meanCNull,{'good','bad'},...
    'Predicted Outcome','True Outcome','Null model',OANull);

% Full
nexttile
show_confusion(meanCFull,{'good','bad'},...
    'Predicted Outcome','True Outcome','Full model',OAFull);

end