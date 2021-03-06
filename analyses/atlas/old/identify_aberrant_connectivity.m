function identify_aberrant_connectivity

%% To do
%{
I SCREWED UP AND DIDN'T RESHIFT SOZ WHEN I DELETED STUFF
- use john's atlas? Need his code to convert mni coordinates to AAL
%}


%% Parameters
which_atlas = 'aal_bernabei_bipolar';
normalize_edges = 0;
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

%% convert sozs to binary matrix
all_regions = atlas_nums';
bin_soz = cell2mat(cellfun(@(x) ismember(all_regions,x),sozs,'uniformoutput',false));
bin_soz = bin_soz';

%% Confirm soz atlas names seem reasonable
if 1
    soz_names = cell(npts,1);
    for i = 1:npts
        %curr_soz = unique(sozs{i});
        %curr_soz_nums = ismember(atlas_nums,curr_soz);
        curr_soz_names = atlas_names(bin_soz(:,i));
        soz_names{i} = curr_soz_names;
    end
    
end

%% Remove cerebellar and not in atlas
%{
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


%ns_ranking_soz(z_score_ns,sozs,atlas_nums,atlas_names,pt_names)
%% Plot all z scores
%{
figure
all_z_score_diff = nan(npts,1);
for i = 1:npts
    curr_soz = bin_soz(:,i);
    z_score_diff = nanmedian(z_score_ns(curr_soz,i))-nanmedian(z_score_ns(:,i));
    all_z_score_diff(i) = z_score_diff;
    %
    plot(i,z_score_diff,'o')
    hold on
    plot(xlim,[0 0],'k--')
    %
    %
    plot(i,nanmedian(z_score_ns(:,i)),'ko')
    hold on
    plot(i,nanmedian(z_score_ns(curr_soz,i)),'ro')
    %
end
sum(all_z_score_diff>0)/sum(~isnan(all_z_score_diff))

%}


%% Plot orders
plot_orders_mats(z_score_ns,bin_soz)



%% Show orders of sozs
%{
% Convert z_score_ns to cell
zscore_C = (num2cell(z_score_ns,1))';

% convert sozs to binary
sozs_bin = cell(npts,1);
nregions = size(z_score_ns,1);
for i = 1:npts
    chs = atlas_nums;
    sozs_bin{i} = ismember(chs,sozs{i});
end

figure
plot_orders(zscore_C,sozs_bin);
%}


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
