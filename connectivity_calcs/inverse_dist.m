function A = inverse_dist(locs)

n = size(locs,1);

A = nan(n,n);
A(logical(eye(n))) = 0;

for i = 1:n
    for j = 1:i-1
        dist = vecnorm(locs(i,:)-locs(j,:));
        inv_dist = 1/dist;
        A(i,j) = inv_dist;
        A(j,i) = inv_dist;
    end
end

% Wrap it for storage
A = wrap_or_unwrap_adjacency_fc_toolbox(A);
       

end