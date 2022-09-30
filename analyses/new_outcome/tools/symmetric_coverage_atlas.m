function [symm_cov_atlas,all_bilateral] = symmetric_coverage_atlas(atlas,locs,lats)

%{
Build atlas where all regions without symmetric coverage are nans
%}

npts = size(atlas,1);

symm_cov_atlas = nan(size(atlas));
dims = ndims(atlas);


% Find the regions with symmetric coverage for each patient
all_bilateral = find_bilateral(atlas,lats,locs);

% confirm that these regions are symmetric
assert(isequal(all_bilateral(strcmp(lats,'L'),:),...
    all_bilateral(strcmp(lats,'R'),:)))

for ip = 1:npts
    curr_bilateral = logical(all_bilateral(:,ip));
    if sum(curr_bilateral) == 0, continue; end
    switch dims
        case 2
            curr_atlas = squeeze(atlas(ip,:));
            curr_atlas(~curr_bilateral) = nan;
            symm_cov_atlas(ip,:) = curr_atlas;
        case 3
            curr_atlas = squeeze(atlas(ip,:,:));
            if size(curr_atlas,1) == size(curr_atlas,2)
                curr_atlas(~curr_bilateral,:) = nan; 
                curr_atlas(:,~curr_bilateral) = nan;
                symm_cov_atlas(ip,:,:) = curr_atlas;
            else
                curr_atlas(~curr_bilateral,:) = nan; 
                symm_cov_atlas(ip,:,:) = curr_atlas;
            end
        case 4
            curr_atlas = squeeze(atlas(ip,:,:,:));
            curr_atlas(~curr_bilateral,:,:) = nan; 
            curr_atlas(:,~curr_bilateral,:) = nan;
            symm_cov_atlas(ip,:,:,:) = curr_atlas;
    end
    

end


end