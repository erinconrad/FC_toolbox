function show_atlas

%% Parameters
which_atlas = 'aal';
gamma = 1; % for community detection

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

%% Remove all nan rows
nan_rows = sum(isnan(avg_atlas),1) == size(avg_atlas,2);
avg_atlas(nan_rows,:) = [];
avg_atlas(:,nan_rows) = [];
names(nan_rows) =[];

%% replace remaining nans with zeros
avg_atlas(isnan(avg_atlas)) = 0;

%% Plot it
figure
set(gcf,'position',[10 10 900 900])
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

%% Detect communities
[Ci,Q]=modularity_und(avg_atlas,gamma);

% Show communities
for i = 1:length(unique(Ci))
    curr_community = Ci == i;
    fprintf('\nCommunity %d:\n',i);
    table(names(curr_community))
end

end