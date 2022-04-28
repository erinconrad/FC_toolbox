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
soz = out.all_soz_bin;
npts = length(soz);
locs = out.all_locs;

%% Calculate default search radius
sr = calculate_default_search_radius(locs); 
% this returns something too small! My intuition is that

%% Compare coverage density for SOZ vs not
if 1
dens = nan(npts,2);
for ip = 1:npts
    
    curr_soz = soz{ip};
    curr_locs = locs{ip};
    
    % get density
    density = estimate_coverage_density(curr_locs,sr);
    % print(gcf,[plot_folder,'density_example'],'-dpng')
    
    dens(ip,:) = [nanmedian(density(curr_soz==1)) nanmedian(density(curr_soz==0))];
    
end

paired_plot(dens,'Density',{'SOZ','non-SOZ'})
end


end