function coeff_stats = sleep_model_bootstrap_stats

rng(0)
nb = 1e3;
ncoeffs = 4;
pt_stats = nan(nb,ncoeffs);
for ib = 1:nb
    
    mout = updated_classifier_may2022([],1,[],[]);
    if ~isfield(mout,'labels'), continue; end
    for ic = 1:ncoeffs
        pt_stats(ib,ic) = mout.model.Coefficients{ic+1,2}; 
    end
end

coeff_stats = nan(ncoeffs,4); % mean, lower CI, higher CI, p
for ic = 1:ncoeffs
    tout = bootstrap_ci_and_p(squeeze(pt_stats(:,ic)));
    coeff_stats(ic,:) = [tout.mean,tout.CI_95,tout.p];
    
end

end