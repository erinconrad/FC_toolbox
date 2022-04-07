function coverage_density_soz

%% File locations
locations = fc_toolbox_locs;
addpath(genpath(locations.script_folder))
results_folder = [locations.main_folder,'results/'];
spike_folder = [results_folder,'analysis/outcome/data/'];
atlas_folder = [results_folder,'analysis/atlas/'];
plot_folder = [results_folder,'analysis/outcome/plots/'];

%% Load out file and get roc stuff
out = load([spike_folder,'main_out.mat']);
out = out.out;

%% Get stuff
rate = out.all_spikes;
soz = out.all_soz_bin;
npts = length(soz);
labels = out.all_labels;
ns = out.all_ns;
fc = out.all_fc;
locs = out.all_locs;

%% Calculate default search radius
sr = calculate_default_search_radius(locs);

%% Compare coverage density for SOZ vs not
if 0
dens = nan(npts,2);
for ip = 1:npts
    
    curr_soz = soz{ip};
    curr_locs = locs{ip};
    
    % get density
    density = estimate_coverage_density(curr_locs,1e3);
    
    dens(ip,:) = [nanmean(density(curr_soz==1)) nanmean(density(curr_soz==0))];
    
end

paired_plot(dens,'Density',{'SOZ','non-SOZ'})
end

%% Classifier to predict 

end