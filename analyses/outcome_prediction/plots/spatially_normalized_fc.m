function spatially_normalized_fc

%% params
which_atlas = 'aal_bernabei'; %'brainnetome';%

%% Locs
locations = fc_toolbox_locs;
script_folder = locations.script_folder;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/outcome/data/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];
atlas_folder = [results_folder,'analysis/atlas/'];

%% Load out file
out = load([out_folder,'main_out.mat']);
out = out.out;

%% get stuff
locs = out.all_locs;
fc = out.all_fc;
soz = out.all_soz_bin;
npts = length(fc);
labels = out.all_labels;
rate = out.all_spikes;

%% Turn soz to logical
soz = cellfun(@logical,soz,'uniformoutput',false);

%% Load atlas file
atlas_out = load([atlas_folder,which_atlas,'.mat']);
atlas_out = atlas_out.out;


%% get atlas stuff
atlas_elec_regions = atlas_out.elecs_atlas;
atlas_nums = atlas_out.atlas_nums;
atlas_elec_labels = atlas_out.elecs_labels;


%% Spatially normalize the FC matrix
resid = fit_distance_model(locs,fc,soz);

%% Get ns of fc and resid
ns = cellfun(@(x) nanmean(x,2),fc,'uniformoutput',false);
ns_resid = cellfun(@(x) nanmean(x,2),resid,'uniformoutput',false);

%% Normalize by anatomy
[z,ana_fc_atlas] = normalize_by_anatomy(ns_resid,atlas_elec_regions,atlas_nums,atlas_elec_labels,soz);



figure
tiledlayout(3,1,'tilespacing','tight','padding','tight')
all = logical(ones(npts,1));

%{
nexttile
do_plot(rate,soz,all,'all','Node strength')
%}

nexttile
do_plot(ns,soz,all,'all','Node strength')

nexttile
do_plot(ns_resid,soz,all,'all','Node strength (distance normalized)')

nexttile
do_plot(z,soz,all,'all','Node strength (distance and anatomically normalized)')


end

function do_plot(rate,soz,which_pts,pt_text,thing_to_plot)

curr_things = rate(which_pts);
curr_soz = soz(which_pts);

% remove nans
old_curr_things = curr_things;
curr_things = cellfun(@(x,y) x(~isnan(x) & ~isnan(y)),curr_things,curr_soz,'uniformoutput',false);
curr_soz = cellfun(@(x,y) y(~isnan(x) & ~isnan(y)),old_curr_things,curr_soz,'uniformoutput',false);

% make sure no nans
for ip = 1:length(curr_things)
    assert(~any(isnan(curr_things{ip})))
    assert(~any(isnan(curr_soz{ip})))
end

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

function [z,atlas] = normalize_by_anatomy(thing,atlas_elec_regions,atlas_nums,atlas_labels,soz)
    
%% Remove ekg from atlas
npts = length(thing);
for ip = 1:npts
    curr_atlas_labels = atlas_labels{ip};
    curr_atlas_elec_regions = atlas_elec_regions{ip};
    ekg = find_non_intracranial(curr_atlas_labels);
    atlas_elec_regions{ip} = curr_atlas_elec_regions(~ekg);
end


nregions = length(atlas_nums);

atlas = nan(nregions,npts);

%% Build atlas
for ip = 1:npts
    curr_thing = thing{ip};
    curr_elecs = atlas_elec_regions{ip};
    curr_soz = soz{ip};
    
    for ir = 1:nregions
        curr_region = atlas_nums(ir);
        
        % get the elecs in that region
        elec_idx =  curr_elecs == curr_region;
        
        % remove SOZ
        if any(elec_idx & curr_soz), continue; end
        
        % get the average thing in those elecs
        atlas(ir,ip) = nanmean(curr_thing(elec_idx));
    end
    
    
end

%% Normalize
z = cell(npts,1);
for ip = 1:npts
    curr_thing = thing{ip};
    curr_elecs = atlas_elec_regions{ip};
    z{ip} = nan(length(curr_thing),1);
    
    for ich = 1:length(curr_thing)
        curr_region = curr_elecs(ich);
        
        if isnan(curr_region), continue; end
        
        curr_atlas_row = atlas_nums == curr_region; 
        z{ip}(ich) = (curr_thing(ich) - nanmean(atlas(curr_atlas_row,:)))./...
            nanstd(atlas(curr_atlas_row,:),[],2);
    end
end

end
