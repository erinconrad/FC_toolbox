function density_soz(all_locs,all_soz)

npts = length(all_locs);

all_dens = [];
all_soz = [];
for ip = 1:npts
    locs = all_locs{ip};
    soz = all_soz{ip};
    
    % make inter-distance matrix
    D = make_interdist_matrix(locs);
    
    
end

end