function out_probs = simple_loo_cv(T,modelspec,distribution)

npts = size(T,1);
out_probs = nan(npts,1);
for ip = 1:npts
    train = [1:ip-1,ip+1:npts];
    test = ip;

    Ttrain = T(train,:);
    Ttest = T(test,:);

    mdl = fitglm(Ttrain,modelspec,'Distribution',distribution);
    prob_test = predict(mdl,Ttest);
    out_probs(ip) = prob_test;

end



end