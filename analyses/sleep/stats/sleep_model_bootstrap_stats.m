function [coeff_names,coeff_stats] = sleep_model_bootstrap_stats(just_gray)


nb = 1e3;
ncoeffs = 4;
pt_stats = nan(nb,ncoeffs);
coeff_names_expected = {'vec_rate_sleep','vec_rate_wake','vec_rate_pre','vec_rate_post'};
for ib = 1:nb
    
    mout = updated_classifier_may2022([],1,[],[],just_gray);
    if ~isfield(mout,'labels'), continue; end
    coeff_names = mout.model.CoefficientNames(2:5);
    assert(isequal(coeff_names,coeff_names_expected))
    
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