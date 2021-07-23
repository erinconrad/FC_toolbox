function ns_stability(pc)

%% Parameters
plotm = 1;
plotf = 1;
spacing = 20;

%% Get file locs
locations = fc_toolbox_locs;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Get network over time
out = net_over_time(pc);

%% Info
ns = out.file(plotf).montage(plotm).ns;
nruns = size(ns,2);
nchs = size(ns,1);
run_center = out.file(plotf).run_center;
clean_labels = out.file(plotf).montage(plotm).labels;


end