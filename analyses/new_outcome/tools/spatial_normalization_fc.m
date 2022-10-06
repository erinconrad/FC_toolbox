function out = spatial_normalization_fc(thing,locs)

%% Make interdistance matrix
D = cellfun(@(x) make_interdist_matrix(x), locs,'UniformOutput',false);

%% Vectorize everything
D_vec = cell2mat(cellfun(@wrap_or_unwrap_adjacency_fc_toolbox,D,'UniformOutput',false));
D_vec(D_vec>1e3) = nan;
thing_vec = cell2mat(cellfun(@wrap_or_unwrap_adjacency_fc_toolbox,thing,'UniformOutput',false));

%% Rat11
bad = isnan(D_vec) | any(isnan(thing_vec),2);
x = D_vec;
y = thing_vec;
x(bad) = [];
y(bad,:) = [];
nfreq = size(y,2);

clear model_out
for jf = 1:nfreq
    model_out(jf).f = fit(x,y(:,jf),'rat11','startpoint',[0.1 0.1 0.1]);
    if 0
        predict_y = (model_out(jf).f.p1 * x + model_out(jf).f.p2)./(x + model_out(jf).f.q1);
        figure
        plot(x,y(:,jf),'o');
        hold on
        plot(x,predict_y,'o');
        pause
    end
end

%% rebuild the matrix of residuals
out = cell(size(thing));
npts = length(thing);
for ip = 1:npts
    
    curr_locs = locs{ip};
    curr_thing = thing{ip};
    
    % make inter-distance matrix
    D = make_interdist_matrix(curr_locs);
    
     % convert to 1D
    D = wrap_or_unwrap_adjacency_fc_toolbox(D);
    curr_thing = wrap_or_unwrap_adjacency_fc_toolbox(curr_thing);
    
    % predict y
    resid = nan(size(curr_thing));
    for jf = 1:nfreq
        predict_y = (model_out(jf).f.p1 * D + model_out(jf).f.p2)./(D + model_out(jf).f.q1);
        
        resid(:,jf) = curr_thing(:,jf) - predict_y;
    end
    
    
    out{ip} = wrap_or_unwrap_adjacency_fc_toolbox(resid);
    
end


end