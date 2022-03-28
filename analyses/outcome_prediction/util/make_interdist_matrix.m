function D = make_interdist_matrix(locs)

nelecs = size(locs,1);
D = nan(nelecs,nelecs);

for i = 1:nelecs
    for j = 1:nelecs-1
        
        locs_i = locs(i,:);
        locs_j = locs(j,:);
        
        d = vecnorm(locs_i-locs_j);
        D(i,j) = d;
        D(j,i) = d;
        
    end
    
end


end