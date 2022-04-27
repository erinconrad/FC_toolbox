function [f,vec_dens,vec_conn] = inner_dens_model(all_density,all_fc,all_soz_bin,all_spikes,max_spikes)

npts = length(all_density);

%% build model
vec_conn = [];
vec_dens = [];
for ip = 1:npts
    conn = (all_fc{ip});
    soz = all_soz_bin{ip};
    spikes = all_spikes{ip};
    dens = all_density{ip};
    
    soz = logical(soz);
    
    % convert to 1D
    dens = wrap_or_unwrap_adjacency_fc_toolbox(dens);
    conn = wrap_or_unwrap_adjacency_fc_toolbox(conn);
    
    % make soz nans
    conn(soz) = nan;
    
    % make those with many spikes nans
    many_spikes = spikes > max_spikes;
    conn(many_spikes) = nan;
    
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