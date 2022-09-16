function all_bilateral = find_bilateral(atlas,lats,locs)

% find regions with bilateral coverage

nregions = size(atlas,2);
npts = size(atlas,1);
dim = ndims(atlas);

all_bilateral = nan(nregions,npts);

% get left and right regions
left_regions = strcmp(lats,'L');
right_regions = strcmp(lats,'R');

for ip = 1:npts
    switch dim
        case 2
            curr_atlas = atlas(ip,:);
        case 3
            curr_atlas = nanmean(atlas(ip,:,:),3);
        case 4
            curr_atlas = nanmean(atlas(ip,:,:,:),[3 4]);
    end
    bilateral_region = zeros(nregions,1);
    
    % Loop over regions
    for ir = 1:nregions
        
        % skip if not a left region
        if left_regions(ir) == 0, continue; end
        
        % find the row of the corresponding right region
        curr_loc = locs{ir};
        right_region = strcmp(locs,curr_loc) & right_regions; % same loc but right
        
        % Skip if there is no corresponding right region
        if sum(right_region) == 0, continue; end % skipping it will leave it zero
        
        % see if both the left and right region have non-empty elements in
        % the atlas
        if ~isnan(curr_atlas(ir)) && ~isnan(curr_atlas(right_region))
            bilateral_region(ir) = 1; % then it has bilateral coverage
            bilateral_region(right_region) = 1;
        end
        %{
        if sum(~isnan(curr_atlas(ir,:))) > 0 && ...
                sum(~isnan(curr_atlas(right_region,:))) > 0
            bilateral_region(ir) = 1; % then it has bilateral coverage
            bilateral_region(right_region) = 1;
        
        end
        %}
    end
    all_bilateral(:,ip) = bilateral_region;
end


end