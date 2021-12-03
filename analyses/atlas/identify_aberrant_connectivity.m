function identify_aberrant_connectivity

%% To do
%{
- use john's atlas? Need his code to convert mni coordinates to AAL
%}


%% Parameters
which_atlas = 'aal';
normalize_edges = 1;
rm_soz = 0;

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
atlas_names = out.atlas_names;
atlas_nums = out.atlas_nums;
sozs = out.sozs;

npts = size(atlas,3);
pt_names = out.pt_names;

%% Remove cerebellar and not in atlas
%
cerebellar = contains(atlas_names,'Cerebelum');
not_in_atlas = strcmp(atlas_names,'NotInAtlas');

atlas = atlas(~cerebellar & ~not_in_atlas,~cerebellar & ~not_in_atlas,:);
atlas_names = atlas_names(~cerebellar & ~not_in_atlas);
atlas_nums = atlas_nums(~cerebellar & ~not_in_atlas);
%}

%% Remove SOZ
% One thing I worry about this is that I might be introducing a bias
% wherein, by deliberately excluding SOZ electrodes from the atlas, they
% are defined to be more abnormal. So I increase my probability of Type I
% errors.
old_atlas = atlas;
if rm_soz
    atlas = remove_soz_from_atlas(atlas,atlas_nums,sozs);
end



%% Convert atlas to z-score
% Note that old_atlas is the one where I deliberately keep SOZ
% but atlas is the one used to generate norms (so I remove them)
zscore = (old_atlas - nanmean(atlas,3))./nanstd(atlas,[],3);

% Get max and min zscores
min_zscore = min(min(min(zscore)));
max_zscore = max(max(max(zscore)));

%% Get node strength
if normalize_edges
    
    % Take the pre-normalized connection matrix and define the node
    % strength to be the sum across the rows here
    z_score_ns = squeeze(nanmean(zscore,2));
else
    
    % Don't pre-normalize. Take the node strength the usual way and then
    % normalize this
    ns = squeeze(nanmean(atlas,1));

    % Z score ns
    z_score_ns = zscore_anything(ns,2);
    
end


ns_ranking_soz(z_score_ns,sozs,atlas_nums,atlas_names,pt_names)


%% Show individual networks
if 0
figure
set(gcf,'position',[10 10 1400 1400])
for ip = 1:npts
    if sum(sum(~isnan(zscore(:,:,ip)))) == 0, continue; end
    turn_nans_gray(zscore(:,:,ip))
    xticks(1:length(atlas_names))
    yticks(1:length(atlas_names))
    
    title(pt_names{ip})
    colorbar
    caxis([min_zscore max_zscore])
    
    % Add the soz locs
    soz = out.sozs{ip};
    soz_indices = (ismember(atlas_nums,soz));
    special_names = add_asterisks_to_chosen(atlas_names,soz_indices);
    xticklabels(special_names)
    yticklabels(special_names)
    
    pause
end
end


end
