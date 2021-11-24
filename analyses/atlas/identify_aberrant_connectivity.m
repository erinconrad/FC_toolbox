function identify_aberrant_connectivity

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
npts = size(atlas,3);

%% Remove cerebellar and not in atlas
cerebellar = contains(names,'Cerebelum');
not_in_atlas = strcmp(names,'NotInAtlas');

atlas = atlas(~cerebellar & ~not_in_atlas,~cerebellar & ~not_in_atlas,:);
names = names(~cerebellar & ~not_in_atlas);


%% Convert atlas to z-score
zscore = (atlas - nanmean(atlas,3))./nanstd(atlas,[],3);

figure
set(gcf,'position',[10 10 1400 1400])
for ip = 1:npts
    if sum(sum(~isnan(zscore(:,:,ip)))) == 0, continue; end
    turn_nans_gray(zscore(:,:,ip))
    xticks(1:length(names))
    yticks(1:length(names))
    xticklabels(names)
    yticklabels(names)
    
    colorbar
    pause
end


end