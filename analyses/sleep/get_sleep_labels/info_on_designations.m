function info_on_designations

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'analysis/sleep/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

% Load out folder
out = load([out_folder,'out.mat']);
out = out.out;

swdes = out.roc_out.swdes;
nawake = 0;
nasleep = 0;

for i = 1:length(swdes)
    nawake = nawake + length(swdes(i).sw.wake);
    nasleep = nasleep + length(swdes(i).sw.sleep);
end

end