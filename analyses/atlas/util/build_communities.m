function build_communities

%% Parameters
atlas = 'aal_bernabei_bipolar';
gamma_search_space = [0.5:0.01:1.5];

%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
atlas_folder = [results_folder,'analysis/atlas/'];
bct_folder= locations.bct;

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
addpath(genpath(bct_folder));

%% Get avg atlas
out = show_atlas(atlas,0);
avg_atlas = out.avg_atlas;
names = out.names;
nregions = length(names);

%% Get canonical netwokrs associated with aal regions
out = assign_aal_to_canon;
% Find indices corresponding to aal regions I have in my dataset
canon_aal_names = out.aal_names;
yeo_names = out.yeo_names;
aal_to_yeo = out.aal_to_yeo;

% remove no canonical network
aal_to_yeo(:,1) = [];

%% Match names
[C,ia,ib] = intersect(names,canon_aal_names);
names = names(ia);
canon_aal_names = canon_aal_names(ib);
assert(isequal(names,canon_aal_names));
aal_to_yeo = aal_to_yeo(ib,:);
avg_atlas = avg_atlas(ia,ia);

%% Detect communities
nsearch = length(gamma_search_space);
canon_assignments = nan(size(avg_atlas,1),size(aal_to_yeo,2),nsearch);
for is = 1:nsearch
    gamma = gamma_search_space(is);
    canon_assignments(:,:,is) = assign_module_to_canon(avg_atlas,gamma,aal_to_yeo);
end

if 1
    imagesc(nanmean(canon_assignments,3))
    yticks(1:nregions)
    yticklabels(names)

end

%% Build community matrix
communities = unique(Ci);
ncommunities = length(communities);
region_community = zeros(nregions,ncommunities);
for i = 1:nregions
    which_comm = Ci(i);
    region_community(i,which_comm) = 1;
end

figure
tiledlayout(1,2)

nexttile
imagesc(region_community)
yticks(1:nregions)
yticklabels(names)

nexttile
imagesc(aal_to_yeo)
yticks(1:nregions)
yticklabels(names)


end