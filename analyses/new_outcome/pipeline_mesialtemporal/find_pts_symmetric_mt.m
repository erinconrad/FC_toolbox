function [which_names,which_pts] = find_pts_symmetric_mt

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
edf_path = [results_folder,'edf_out/'];
data_folder = [locations.main_folder,'data/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

which_names = {};
which_pts = [];

%% Load pt folder
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

for i = 1:length(pt)
    name = pt(i).name;
    if exist([edf_path,name,'/summ.mat'],'file') == 0
        continue
    end
    summ = load([edf_path,name,'/summ.mat']);
    labels = summ.out.labels;
    if sum(contains(labels,'L'))>0 && sum(contains(labels,'R'))>0
        which_pts = [which_pts;i];
        which_names = [which_names;name];
    end

end

T = table(which_pts,which_names);
writetable(T,[edf_path,'symmetric_mt_patients.csv']);

end