function [f,vec_dens,vec_conn] = inner_dens_model(all_density,all_fc)

npts = length(all_density);

%% build model
vec_conn = [];
vec_dens = [];
for ip = 1:npts
    conn = (all_fc{ip});
    dens = all_density{ip};
    
    
    % convert to 1D
    dens = wrap_or_unwrap_adjacency_fc_toolbox(dens);
    conn = wrap_or_unwrap_adjacency_fc_toolbox(conn);
    
    % add to vector
    vec_dens = [vec_dens;dens];
    vec_conn = [vec_conn;conn];

end

%% Remove nans
bad = isnan(vec_dens) | isnan(vec_conn);

% do a model
x = vec_dens;
y = vec_conn;
x(bad) = [];
y(bad) = [];
f = fitlm(x,y);