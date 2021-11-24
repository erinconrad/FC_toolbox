function show_atlas

%% Parameters
which_atlas = 'aal';

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];
bct_folder= locations.bct;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Load atlas
out = load([atlas_folder,which_atlas,'.mat']);
out = out.out;

atlas = out.atlas;
names = out.atlas_names;

%% Remove cerebellar and not in atlas
cerebellar = contains(names,'Cerebelum');
not_in_atlas = strcmp(names,'NotInAtlas');

atlas = atlas(~cerebellar & ~not_in_atlas,~cerebellar & ~not_in_atlas,:);
names = names(~cerebellar & ~not_in_atlas);
avg_atlas = nanmean(atlas,3);

%% Plot it
turn_nans_gray(avg_atlas)
xticks(1:length(names))
yticks(1:length(names))
xticklabels(names)
yticklabels(names)
colorbar

%% Node strength, sort
ns = nansum(avg_atlas,2);
[ns,I] = sort(ns,'descend');
table(ns,names(I))

end