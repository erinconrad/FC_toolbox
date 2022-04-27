function maxD = max_inter_elec_dist(all_locs)

npts = length(all_locs);
all_max_D = nan(npts,1);

for ip = 1:npts
    locs = all_locs{ip};
    % Fix for weird locs
    locs(any(locs > 1e5,2),:) = nan;
    D = make_interdist_matrix(locs);
    
    maxD = max(D(:));
    all_max_D(ip) = maxD;
end

maxD = max(all_max_D);

end