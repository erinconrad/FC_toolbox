function ns_stability(pc)

%% Parameters
plotm = 1;
spacing = 20;

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Get network over time
out = net_over_time(pc);
out = reconcile_files(out);


%% Info
ns = out.montage(plotm).ns;
nruns = size(ns,2);
nchs = size(ns,1);
clean_labels = out.montage(plotm).labels;

%% Get ns stability
% Mean ns
mean_ns = nanmean(ns,2);
nst = corr(mean_ns,ns,'Type','Spearman','rows','pairwise');

%% Plot
plot(nst);


end