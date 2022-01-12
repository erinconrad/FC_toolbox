function A = make_dist_network(locs)

nchs = size(locs,1);

A = nan(nchs,nchs);
for i = 1:nchs
    for j = 1:i-1
        dist = vecnorm(locs(i,:)-locs(j,:));
        metric = 1/dist^2;
        A(i,j) = metric;
        A(j,i) = metric;
    end
end

end