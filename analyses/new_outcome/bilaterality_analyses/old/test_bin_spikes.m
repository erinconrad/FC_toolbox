function test_bin_spikes(T)

temporal = strcmp(T.soz_locs,'temporal');
T(~temporal,:) = [];

newT = table(T.soz_lats,T.('spikes bipolar sleep'));
newT.response = true(size(newT.Var1,1),1);
newT.response(strcmp(newT.Var1,'right')) = false;
newT.predictor = newT.Var2;
newT.bin_predictor = ones(size(newT.predictor,1),1);
newT.bin_predictor(newT.predictor<0) = 0;
newT.logical_predictor = true(size(newT.predictor,1),1);
newT.logical_predictor(newT.bin_predictor == 0) = false;

mdl1 = fitglm(newT,'response ~ predictor','Distribution','binomial')
mdl2 = fitglm(newT,'response ~ bin_predictor','Distribution','binomial')
mdl3 = fitglm(newT,'response ~ logical_predictor','Distribution','binomial')

npts = size(newT,1);
scores = nan(npts,1);
for ip = 1:npts
    train = [1:ip-1,ip+1:npts];
    mdl = fitglm(newT(train,:),'response ~ logical_predictor','Distribution','binomial');
    scores(ip) = predict(mdl,newT(ip,:));
end

[~,~,~,auc] = perfcurve(newT.response,scores,'true')

end