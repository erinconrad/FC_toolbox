function ranking_by_localization

%% Parameters
max_spikes = 1/3600; % max spikes/elecs/min to include in model
thing_to_plot = 'ns_resid';
which_atlas = 'aal_bernabei'; %'brainnetome';%

rl_min_spikes = 0;
plot_type = 'scatter';
nblocks = 6;
myColours = [0.1660, 0.540, 0.1880;...
0.4940, 0.1840, 0.5560;...    
0.8500, 0.4250, 0.0980;...
    0.9290 0.6940 0.1250];



locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
atlas_folder = [results_folder,'analysis/atlas/'];


%% Load out file and get roc stuff
out = load([out_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
rl = out.all_rl;
ns = out.all_ns;
soz = out.all_soz_bin;
soz_loc = out.all_soz_locs;
npts = length(soz);
labels = out.all_labels;
locs = out.all_locs;
fc = out.all_fc;

%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

%% get atlas stuff
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;
spikes_atlas = atlas_out.spikes_atlas;
atlas_nums = atlas_out.atlas_nums;
atlas_names = atlas_out.atlas_names;
ns_atlas = atlas_out.atlas;
ns_atlas = squeeze(nanmean(ns_atlas,2));

%% Normalize spike rates and ns
z_rates = cell(npts,1);
z_ns = cell(npts,1);
for ip = 1:npts
    z_rates{ip} = normalize_spike_rates(labels{ip},atlas_elec_labels{ip},...
        atlas_elec_regions{ip},spikes_atlas,rate{ip},atlas_nums);
    
    z_ns{ip} = normalize_spike_rates(labels{ip},atlas_elec_labels{ip},...
        atlas_elec_regions{ip},ns_atlas,ns{ip},atlas_nums);
end


%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);

%% For RL, set as nans those without enough spikes
rl = cellfun(@(x,y) make_non_spikey_nan(x,y,rl_min_spikes), rl, rate,'uniformoutput',false);

%% Spatially normalize the FC matrix
resid = fit_distance_model(locs,fc,soz,rate,max_spikes,plot_folder);
ns_resid = cellfun(@(x) nanmean(x,2),resid,'uniformoutput',false);

%% Separate patients by localization
mf = contains(soz_loc,'multifocal') | contains(soz_loc,'diffuse');
mt = strcmp(soz_loc,'mesial temporal');
tn = strcmp(soz_loc,'temporal neocortical');
oc = strcmp(soz_loc,'other cortex');

%% Decide what to plot
switch thing_to_plot
    case 'ns'
        thing = ns;
    case 'rate'
        thing = rate;
    case 'rl'
        thing = rl;
    case 'rate_norm'
        thing = z_rates;
    case 'ns_norm'
        thing = z_ns;
    case 'ns_resid'
        thing = ns_resid;
        
end

figure
set(gcf,'position',[10 10 1100 600])
tiledlayout(3,2,'tilespacing','tight','padding','tight')

%% SOZ spike rate ranking
nexttile([1 2])
all = logical(ones(npts,1));
do_plot(thing,soz,all,'all',thing_to_plot)

nexttile
do_plot(thing,soz,mt,'mesial temporal',thing_to_plot)

nexttile
do_plot(thing,soz,tn,'temporal neocortical',thing_to_plot)

nexttile
do_plot(thing,soz,oc,'other cortex',thing_to_plot)

nexttile
do_plot(thing,soz,mf,'multifocal',thing_to_plot)

print(gcf,[plot_folder,thing_to_plot,'_rankings'],'-dpng')

end

function do_plot(rate,soz,which_pts,pt_text,thing_to_plot)

curr_things = rate(which_pts);
curr_soz = soz(which_pts);

% remove nans
old_curr_things = curr_things;
curr_things = cellfun(@(x,y) x(~isnan(x) & ~isnan(y)),curr_things,curr_soz,'uniformoutput',false);
curr_soz = cellfun(@(x,y) y(~isnan(x) & ~isnan(y)),old_curr_things,curr_soz,'uniformoutput',false);

stats_out = outcome_plot_orders(curr_things,curr_soz);
hold on
xticklabels([])
xlabel('Patient')
ylabel(sprintf('Electrode %s rank',thing_to_plot))
set(gca,'fontsize',15)
title(sprintf('SOZ %s ranking for %s onsets',thing_to_plot,pt_text))
xl = xlim;
yl=ylim;
text(mean(xl),yl(2),sprintf('median rank = %1.1f',stats_out.median_rank),...
    'horizontalalignment','center','verticalalignment','top','fontsize',15)

end

function x = make_non_spikey_nan(x,y,rl_min_spikes)

x(y<rl_min_spikes) = nan;

end

function z = normalize_spike_rates(sp_labels,atlas_labels,regions,atlas,spikes,atlas_nums)
    
%% Remove ekg from atlas
ekg = find_non_intracranial(atlas_labels);
atlas_labels(ekg) = [];
regions(ekg) = [];

%% Reconcile labels
assert(isequal(sp_labels,atlas_labels)) 

z = nan(length(spikes),1);

% Loop over electrodes
for ich = 1:length(spikes)
    
    % get the region of this electrode
    curr_region = regions(ich);
    
    % get the index of this
    curr_idx = find(atlas_nums == curr_region);
    
    if isempty(curr_idx), z(ich) = nan; continue; end
    
    % normalize by atlas (across all patients for that region)
    %z(ich) = (spikes(ich)-nanmedian(atlas(curr_idx,:)))./iqr(atlas(curr_idx,:));
    z(ich) = (spikes(ich)-nanmean(atlas(curr_idx,:)))./nanstd(atlas(curr_idx,:));
end

if 0
    table(sp_labels,spikes,z)
end

end
