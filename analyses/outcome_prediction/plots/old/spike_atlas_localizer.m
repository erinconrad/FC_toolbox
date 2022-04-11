function spike_atlas_localizer

%% Parameters
which_atlas = 'aal_bernabei';%'brainnetome';%'aal_bernabei';

%% File locations
locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
spike_folder = [results_folder,'analysis/outcome/data/'];
atlas_folder = [results_folder,'analysis/atlas/'];


%% Load out file and get roc stuff
out = load([spike_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
rl = out.all_rl;
soz = out.all_soz_bin;
loc = out.all_soz_locs;
npts = length(soz);
labels = out.all_labels;

%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);

%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;

%% get atlas stuff
atlas_elec_labels = atlas_out.elecs_labels;
atlas_elec_regions = atlas_out.elecs_atlas;
spikes_atlas = atlas_out.spikes_atlas;
atlas_nums = atlas_out.atlas_nums;

%% Normalize spike rates
z = cell(npts,1);
for ip = 1:npts
    z{ip} = normalize_spike_rates(labels{ip},atlas_elec_labels{ip},...
        atlas_elec_regions{ip},spikes_atlas,rate{ip},atlas_nums);
end

%% Plot

figure
tiledlayout(2,1)

%% SOZ spike rate ranking
nexttile
do_plot(rate,soz)

%% SOZ spike rate ranking
nexttile
do_plot(z,soz)


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


function do_plot(rate,soz)

stats_out = plot_orders(rate,soz);
hold on
xticklabels([])
xlabel('Patient')
ylabel('Electrode spike rate rank')
set(gca,'fontsize',15)
title(sprintf('SOZ spike rate ranking'))
xl = xlim;
yl=ylim;
text(mean(xl),yl(2),sprintf('median rank = %1.1f',stats_out.median_rank),...
    'horizontalalignment','center','verticalalignment','top','fontsize',15)

end