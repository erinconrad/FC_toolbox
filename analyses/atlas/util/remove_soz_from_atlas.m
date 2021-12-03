function atlas = remove_soz_from_atlas(atlas,atlas_nums,sozs)

old_atlas = atlas;

for ip = 1:size(atlas,3)
    curr_atlas = atlas(:,:,ip);
    
    % get soz atlas nums
    curr_soz = sozs{ip};
    
    % remove nans
    curr_soz(isnan(curr_soz)) = [];
    
    % get the corresponding indices of the atlas matrix
    [lia,soz_indices] = ismember(curr_soz,atlas_nums);
    
    if any(lia==0), error('cannot find soz atlas number in atlas nums'); end
    
    % get the unique soz indices
    soz_indices = unique(soz_indices);
    
    % make all connections to the soz nans (because I
    % suspect they are abnormal)
    curr_atlas(soz_indices,:) = nan;
    curr_atlas(:,soz_indices) = nan;
    
    atlas(:,:,ip) = curr_atlas;
    
end

if 0
    % compare old to new atlas after removing soz (end up looking nearly
    % identical)
    figure
    tiledlayout(1,3)
    nexttile
    turn_nans_gray(nanmean(old_atlas,3))
    title('Old')
    
    nexttile
    turn_nans_gray(nanmean(atlas,3))
    title('New')
    
    nexttile
    plot(old_atlas(:),atlas(:),'o')
    r = corr(old_atlas(:),atlas(:),'rows','pairwise');
    title(sprintf('Correlation between old and new r = %1.2f',r));
end



end