function [symm_cov_atlas,all_bilateral] = build_symmetric_coverage_atlas(atlas,locs,lats)

npts = size(atlas,3);

symm_cov_atlas = nan(size(atlas));

% Find the regions with symmetric coverage for each patient
all_bilateral = find_symmetric_coverage(atlas,lats,locs);

% confirm that these regions are symmetric
assert(isequal(all_bilateral(strcmp(lats,'L'),:),...
    all_bilateral(strcmp(lats,'R'),:)))

for ip = 1:npts
    curr_bilateral = logical(all_bilateral(:,ip));
    if sum(curr_bilateral) == 0, continue; end
    curr_atlas = atlas(:,:,ip);
    curr_atlas(~curr_bilateral,:) = nan;
    curr_atlas(:,~curr_bilateral) = nan;
    symm_cov_atlas(:,:,ip) = curr_atlas;

end


end