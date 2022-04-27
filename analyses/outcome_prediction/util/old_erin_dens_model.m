function out = erin_dens_model(nout,max_spikes,plot_folder,sr)
unpack_any_struct(nout)


%% Nodal level:
%{
% Get all electrode densities
% default search radius
sr = calculate_default_search_radius(all_locs);

% get all densities
all_dens = cell(npts,1);
for i = 1:npts
    locs = all_locs{i};
    density = estimate_coverage_density(locs,sr);
    all_dens{i} = density;
end

% Concatenate all densities and ns into vectors
vec_dens = [];
vec_ns = [];
vec_ns_sq = [];
for i = 1:npts
    vec_dens = [vec_dens;all_dens{i}];
    
    % get ns
    fc = (all_fc{i});
    fc_sq = (all_fc{i}).^2;
    ns = nansum(fc,2);
    ns_sq = nansum(fc_sq,2);
    vec_ns = [vec_ns;ns];
    vec_ns_sq = [vec_ns_sq;ns_sq];
end

% Plot
if 0
figure
nexttile
plot(vec_dens,vec_ns,'o')
xlabel('Density')
ylabel('NS (r)')

nexttile
plot(vec_dens,vec_ns_sq,'o')
xlabel('Density')
ylabel('NS (r^2)')
end
%}

%% Edge level
% For each electrode i, estimate the density from every other electrode j

%all_fc = cellfun(@(x) x.^2, all_fc,'uniformoutput',false);

npts = length(all_locs);
% default search radius

if ~exist('sr','var') || isempty(sr)
    sr = calculate_default_search_radius(all_locs);
    %sr = max_inter_elec_dist(all_locs);
end

vec_conn = [];
vec_dens = [];
which_pt = [];

for ip = 1:npts
    locs = all_locs{ip};
    conn = (all_fc{ip});
    soz = all_soz_bin{ip};
    spikes = all_spikes{ip};
    
    soz = logical(soz);
    
    % Fix for weird locs
    locs(any(locs > 1e5,2),:) = nan;
    
    % make inter-distance matrix
    D = make_interdist_matrix(locs);
    
    % Make density matrix
    dens = interdistance_to_density_matrix(D,sr);
    
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
    which_pt = [which_pt;repmat(ip,length(dens),1)];
    
end

%% Remove nans and very large distances
bad = isnan(vec_dens) | isnan(vec_conn);

% do a model
x = vec_dens;
y = vec_conn;
x(bad) = [];
y(bad) = [];
f = fitlm(x,y);
est = f.Coefficients.Estimate;
g.p1 = est(1);
g.p2 = est(2);

% rebuild the matrix of residuals
out = cell(size(all_fc));
for ip = 1:npts
    
    locs = all_locs{ip};
    fc = all_fc{ip};
    
    % make inter-distance matrix
    D = make_interdist_matrix(locs);
    
    % Make density matrix
    dens = interdistance_to_density_matrix(D,sr);
    
     % convert to 1D
    dens = wrap_or_unwrap_adjacency_fc_toolbox(dens);
    fc = wrap_or_unwrap_adjacency_fc_toolbox(fc);
    
    % predict y
    predict_y = g.p1 + dens*g.p2;
    
    resid = fc - predict_y;
    
    
    out{ip} = wrap_or_unwrap_adjacency_fc_toolbox(resid);
    
end

if exist('plot_folder','var') && ~isempty(plot_folder)
    ip = 90;
    
    
    figure
    set(gcf,'position',[10 10 1000 650])
    tiledlayout(2,2,'tilespacing','tight','padding','tight')
    
    nexttile
    conn = all_fc{ip};
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
    dens = interdistance_to_density_matrix(D,sr);
    turn_nans_gray(dens)
    c = colorbar;
    ylabel(c,'Density (mm)')
    set(gca,'fontsize',15)
    xticklabels([])
    yticklabels([])
    xlabel('Electrode')
    ylabel('Electrode')
    
    nexttile
    plot(vec_dens,vec_conn,'o')
    hold on
    temp_x = [min(x):0.5:max(x)];
    %temp_y = (f.p1 * temp_x + f.p2)./(temp_x + f.q1);
    temp_y = (g.p1 + temp_x*g.p2);
    plot(temp_x,temp_y,'k','linewidth',3)
    ylim([-1 1])
    xlim([min(x) max(x)])
    xlabel('Density (mm)')
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
    
    print(gcf,[plot_folder,'dens_model'],'-dpng')
end


end