function [coeff_names,coeff_stats] = sleep_model_bootstrap_stats(just_gray,pre_implant)


nb = 1e3; % number of bootstrap samples, should be 1e3

%% Establish the fixed effects
ncoeffs = 4;
pt_stats = nan(nb,ncoeffs);
coeff_names_expected = {'vec_rate_sleep','vec_rate_wake','vec_rate_pre','vec_rate_post'};

% Loop over bootstrap samples
for ib = 1:nb

    % wrap in a while loop to retry if funny errors
    while 1
        %% Do the classifier
        mout = updated_classifier_may2022([],1,[],[],just_gray,pre_implant);
        %{ 
        first argument [] means don't leave any patients out (do bootstrap
        sampling instead. Second argument 1 means do mixed effects model. 3rd
        argument [] means not doing the wake vs sleep analysis. 4th argment []
        means full duration. just_gray indicates whether to do all electrodes
        or only gray matter.

        %}
    
        if ~(isfield(mout,'labels'))
        % If this happens, it means model failed. I don't like this! Retry
            continue
        else
            % accept this model
            break
        end
    end
    
    coeff_names = mout.model.CoefficientNames(2:5);
    assert(isequal(coeff_names,coeff_names_expected))
    
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