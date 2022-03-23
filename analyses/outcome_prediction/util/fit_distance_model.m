function out = fit_distance_model(all_locs,all_conn,all_soz)

npts = length(all_locs);
vec_dist = [];
vec_conn = [];
which_pt = [];
for ip = 1:npts
    locs = all_locs{ip};
    conn = all_conn{ip};
    soz = all_soz{ip};
    
    % make inter-distance matrix
    D = make_interdist_matrix(locs);
    
    % convert to 1D
    D = wrap_or_unwrap_adjacency_fc_toolbox(D);
    conn = wrap_or_unwrap_adjacency_fc_toolbox(conn);
    
    % make soz nans
    conn(soz) = nan;
    
    % add to vector
    vec_dist = [vec_dist;D];
    vec_conn = [vec_conn;conn];
    which_pt = [which_pt;repmat(ip,length(D),1)];
    
end

%% Remove nans and very large distances
bad = isnan(vec_dist) | isnan(vec_conn) | vec_dist > 1e3;
vec_dist(bad) = [];
vec_conn(bad) = [];

% show it
if 0
    figure
    plot(vec_dist,vec_conn,'o')
end

% do a model
x = vec_dist;
y = vec_conn;
f = fit(x,y,'rat11','startpoint',[0.1 0.1 0.1]);

% rebuild predicted y
predict_y = (f.p1 * x + f.p2)./(x + f.q1);

% calculate the residuals
resid = y-predict_y;

% rebuild the matrix
all_resid = nan(length(bad),1);
all_resid(~bad) = resid;

out = cell(size(all_conn));
for ip = 1:npts
    idx = which_pt == ip;
    
    curr_resid = all_resid(idx);
    out{ip} = wrap_or_unwrap_adjacency_fc_toolbox(curr_resid);
    
end


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