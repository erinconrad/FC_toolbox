function [coeff_names,coeff_stats] = sleep_model_bootstrap_stats(just_gray)


nb = 1e3; % number of bootstrap samples, should be 1e3

%% Establish the fixed effects
coeff_names_expected = {'vec_rate_sleep','vec_rate_wake','vec_rate_pre',...
    'vec_rate_post','vec_mri_lesional','vec_concordant_loc','vec_concordant_lat'};
ncoeffs = length(coeff_names_expected);
pt_stats = nan(nb,ncoeffs);

% Loop over bootstrap samples
for ib = 1:nb

    % wrap in a while loop to retry if funny errors
    while 1
        %% Do the classifier
        mout = classifier_with_preimplant([],[],[],just_gray,0);
        %{ 
        first argument [] means don't leave any patients out (do bootstrap
        sampling instead. 2nd
        argument [] means not doing the wake vs sleep analysis. 3rd argment []
        means full duration. just_gray indicates whether to do all electrodes
        or only gray matter.

        %}
    
        if ~(isfield(mout,'labels'))
        % If this happens, it means model failed. Retry
            continue
        else
            assert(~isnan(mout.AUC)) % I don't want any nans because they will mess with my bootstrap analysis
            % accept this model
            break
        end
    end
    
    coeff_names = mout.model.CoefficientNames(2:8);
   % assert(isequal(coeff_names,coeff_names_expected)) % dumb error
    
    % Get individual model coefficients
    for ic = 1:ncoeffs
        pt_stats(ib,ic) = mout.model.Coefficients{ic+1,2}; 
    end
end

%% Bootstap CI and p-values from individual model coefficients
coeff_stats = nan(ncoeffs,4); % mean, lower CI, higher CI, p
for ic = 1:ncoeffs
    tout = bootstrap_ci_and_p(squeeze(pt_stats(:,ic)));
    coeff_stats(ic,:) = [tout.mean,tout.CI_95,tout.p];
    
end

end