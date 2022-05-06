function out = erin_dens_model(nout,max_spikes,plot_folder)

%% Params
poss_sr = [5:5:100];
nsr = length(poss_sr);

%% Unpack struct
all_locs = nout.all_locs;
all_fc = nout.all_fc;
all_soz_bin = nout.all_soz_bin;
all_spikes = nout.all_spikes;
npts = length(all_locs);


%% fix weird locs
for ip = 1:npts
    locs = all_locs{ip};
    locs(any(locs > 1e5,2),:) = nan;
    all_locs{ip} = locs;
end

%% Find optimal sr
all_r2 = nan(nsr,1);
best_sr = nan;
best_f = nan;
best_vec_dens = nan;
best_vec_conn = nan;
for isr = 1:nsr
    sr = poss_sr(isr);
    [f,vec_dens,vec_conn] = dens_model_specific_sr(all_locs,all_fc,all_soz_bin,all_spikes,max_spikes,sr);
    all_r2(isr) = f.Rsquared.Ordinary;
    
    if all_r2(isr) == max(all_r2)
        best_f = f;
        best_sr = sr;
        best_vec_dens = vec_dens;
        best_vec_conn = vec_conn;
    end
end

% Use the best model (highest R2)
f = best_f;
sr = best_sr;
vec_dens = best_vec_dens;
vec_conn = best_vec_conn;
est = f.Coefficients.Estimate;
g.p1 = est(1);
g.p2 = est(2);

% rebuild the matrix of residuals
resid = cell(size(all_fc));
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
    
    tresid = fc - predict_y;
 
    resid{ip} = wrap_or_unwrap_adjacency_fc_toolbox(tresid);
    
end

out.resid = resid;
out.f = f;
out.sr = sr;
out.vec_dens = vec_dens;
out.vec_conn = vec_conn;
out.all_fc = all_fc;
out.all_locs = all_locs;
out.g = g;
out.poss_sr = poss_sr;
out.all_r2 = all_r2;

if exist('plot_folder','var') && ~isempty(plot_folder)
    
    save([plot_folder,'dens_model.mat'],'out');
end


end