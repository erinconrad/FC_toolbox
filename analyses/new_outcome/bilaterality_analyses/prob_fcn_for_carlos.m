function prob = prob_fcn_for_carlos(L,R,params)

pcaWeights = params.pcaWeights;
pcaCoefficients = params.pcaCoefficients;
pcaCenters = params.pcaCenters;
coef = params.coef;

% Calculate AI
AI = (L-R)/(L+R);

% Do PCA (do I need this? Try to remove it)
x = (AI-pcaCenters).*pcaWeights*pcaCoefficients;

% Do logistic function
prob = 1/(1+exp(-coef(1)+coef(2)*x));

end