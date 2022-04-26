function [out,f] = fit_distance_model(all_locs,all_conn,all_soz,all_spikes,max_spikes,plot_folder)

npts = length(all_locs);
vec_dist = [];
vec_conn = [];
which_pt = [];
for ip = 1:npts
    locs = all_locs{ip};
    conn = (all_conn{ip});
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

% do a model
x = vec_dist;
y = vec_conn;
x(bad) = [];
y(bad) = [];
f = fit(x,y,'rat11','startpoint',[0.1 0.1 0.1]);

% rebuild the matrix of residuals
out = cell(size(all_conn));
for ip = 1:npts
    
    locs = all_locs{ip};
    conn = all_conn{ip};
    
    % make inter-distance matrix
    D = make_interdist_matrix(locs);
    
     % convert to 1D
    D = wrap_or_unwrap_adjacency_fc_toolbox(D);
    conn = wrap_or_unwrap_adjacency_fc_toolbox(conn);
    
    % predict y
    predict_y = (f.p1 * D + f.p2)./(D + f.q1);
    
    resid = conn - predict_y;
    
    
    out{ip} = wrap_or_unwrap_adjacency_fc_toolbox(resid);
    
end


% show it
if 0
    ip = 60;
    
    
    figure
    set(gcf,'position',[10 10 1000 650])
    tiledlayout(2,2,'tilespacing','tight','padding','tight')
    
    nexttile
    conn = all_conn{ip};
    turn_nans_gray(conn)
    c = colorbar;
    caxis([-1 1])
    ylabel(c,'Pearson correlation')
    set(gca,'fontsize',15)
    xticklabels([])
    yticklabels([])
    xlabel('Electrode')
    ylabel('Electrode')
    
    nexttile
    D = make_interdist_matrix(all_locs{ip});
    turn_nans_gray(D)
    c = colorbar;
    ylabel(c,'Distance (mm)')
    set(gca,'fontsize',15)
    xticklabels([])
    yticklabels([])
    xlabel('Electrode')
    ylabel('Electrode')
    
    nexttile
    plot(vec_dist,vec_conn,'o')
    hold on
    temp_x = [min(x):0.5:max(x)];
    temp_y = (f.p1 * temp_x + f.p2)./(temp_x + f.q1);
    plot(temp_x,temp_y,'k','linewidth',3)
    ylim([-1 1])
    xlim([min(x) max(x)])
    xlabel('Distance (mm)')
    ylabel('Correlation (r)')
    set(gca,'fontsize',15)
    
    nexttile
    conn = out{ip};
    turn_nans_gray(conn)
    c = colorbar;
    caxis([-1 1])
    ylabel(c,'Normalized correlation (model residuals)')
    set(gca,'fontsize',15)
    xticklabels([])
    yticklabels([])
    xlabel('Electrode')
    ylabel('Electrode')
    
    print(gcf,[plot_folder,'dist_model'],'-dpng')
end


end


