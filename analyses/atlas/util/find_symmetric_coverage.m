function all_bilateral = find_symmetric_coverage(atlas,lats,locs)

nregions = size(atlas,1);
npts = size(atlas,3);

all_bilateral = nan(nregions,npts);

% get left and right regions
left_regions = strcmp(lats,'L');
right_regions = strcmp(lats,'R');

for ip = 1:npts
    curr_atlas = atlas(:,:,ip);
    bilateral_region = zeros(nregions,1);
    
    % Loop over regions
    for ir = 1:nregions
        
        % skip if not a left region
        if left_regions(ir) == 0, continue; end
        
        % find the row of the corresponding right region
        curr_loc = locs{ir};
        right_region = strcmp(locs,curr_loc) & right_regions;
        
        % Skip if there is no corresponding right region
        if sum(right_region) == 0, continue; end
        
        % see if both the left and right region have non-empty elements in
        % the atlas
        if sum(~isnan(curr_atlas(ir,:))) > 0 && ...
                sum(~isnan(curr_atlas(right_region,:))) > 0
            bilateral_region(ir) = 1;
            bilateral_region(right_region) = 1;
        
        end
    end
    all_bilateral(:,ip) = bilateral_region;
end


end