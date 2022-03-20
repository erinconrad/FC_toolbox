function fit_distance_model(all_locs,all_conn)

npts = length(all_locs);
vec_dist = [];
vec_conn = [];
for ip = 1:npts
    locs = all_locs{ip};
    conn = all_conn{ip};
    
    % make inter-distance matrix
    D = make_interdist_matrix(locs);
    
    % add to vector
    vec_dist = [vec_dist;D(:)];
    vec_conn = [vec_conn;conn(:)];
    
end

% show it
if 1
    figure
    plot(vec_dist,vec_conn,'o')
end

% do a model
f = fit(vec_dist,vec_conn,'rat11');

end


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