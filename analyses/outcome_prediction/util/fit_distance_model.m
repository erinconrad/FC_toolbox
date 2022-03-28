function [out,f] = fit_distance_model(all_locs,all_conn,all_soz,all_spikes,max_spikes,plot_folder)

npts = length(all_locs);
vec_dist = [];
vec_conn = [];
which_pt = [];
for ip = 1:npts
    locs = all_locs{ip};
    conn = all_conn{ip};
    soz = all_soz{ip};
    spikes = all_spikes{ip};
    
    % make inter-distance matrix
    D = make_interdist_matrix(locs);
    
    % convert to 1D
    D = wrap_or_unwrap_adjacency_fc_toolbox(D);
    conn = wrap_or_unwrap_adjacency_fc_toolbox(conn);
    
    % make soz nans
    conn(soz) = nan;
    
    % make those with many spikes nans
    many_spikes = spikes > max_spikes;
    conn(many_spikes) = nan;
    
    % add to vector
    vec_dist = [vec_dist;D];
    vec_conn = [vec_conn;conn];
    which_pt = [which_pt;repmat(ip,length(D),1)];
    
end

%% Remove nans and very large distances
bad = isnan(vec_dist) | isnan(vec_conn) | vec_dist > 1e3;
vec_dist(bad) = [];
vec_conn(bad) = [];



% do a model
x = vec_dist;
y = vec_conn;
f = fit(x,y,'rat11','startpoint',[0.1 0.1 0.1]);

% rebuild predicted y
predict_y = (f.p1 * x + f.p2)./(x + f.q1);

% show it
if 0
    figure
    plot(vec_dist,vec_conn,'o')
    hold on
    temp_x = [min(x):0.5:max(x)];
    temp_y = (f.p1 * temp_x + f.p2)./(temp_x + f.q1);
    plot(temp_x,temp_y,'linewidth',3)
    ylim([min(y) max(y)])
    xlim([min(x) max(x)])
    xlabel('Distance (mm)')
    ylabel('Correlation (r^2)')
    set(gca,'fontsize',15)
    print(gcf,[plot_folder,'dist_model'],'-dpng')
end

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


