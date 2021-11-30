function identify_aberrant_connectivity

%% To do
%{
- Think about whether to z-score before or after doing node strength
- remove soz elecs from atlas
- use john's atlas? Need his code to convert mni coordinates to AAL
%}


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
atlas_names = out.atlas_names;
atlas_nums = out.atlas_nums;

npts = size(atlas,3);
pt_names = out.pt_names;

%% Remove cerebellar and not in atlas
%{
cerebellar = contains(names,'Cerebelum');
not_in_atlas = strcmp(names,'NotInAtlas');

atlas = atlas(~cerebellar & ~not_in_atlas,~cerebellar & ~not_in_atlas,:);
atlas_names = atlas_names(~cerebellar & ~not_in_atlas);
%}


%% Convert atlas to z-score
zscore = (atlas - nanmean(atlas,3))./nanstd(atlas,[],3);

%% Get ns from this
z_score_ns = squeeze(nanmean(zscore,2));

%% Get ns of soz
if 1
    figure
    set(gcf,'position',[10 10 1400 1400])
    for ip = 1:npts
        if sum(sum(~isnan(zscore(:,:,ip)))) == 0, continue; end
        nums = 1:length(atlas_names);
        plot(nums,z_score_ns(:,ip),'o')
        hold on
        
        % Add the soz locs
        soz = out.sozs{ip};
        soz_indices = (ismember(atlas_nums,soz));
        
        plot(nums(soz_indices),z_score_ns(soz_indices,ip),'ro');
        title(pt_names{ip})
        pause
        hold off
    end

end


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
