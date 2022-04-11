function W = hanley_mcneil(models,true)

nmodels = length(models);
z = nan(nmodels,nmodels);
p = nan(nmodels,nmodels);
W = nan(nmodels,1);
SE = nan(nmodels,1);
r = nan(nmodels,nmodels);

for im = 1:nmodels
    classi = models{im};
    
    % Get AUC and standard error of model i
    [Wi,SEi] = manual_roc(true,classi);
    W(im) = Wi;
    SE(im) = SEi;
    
    for jm = 1:im-1
        
        classj = models{jm};
    
        % Get AUC and standard error of model i
        [Wj,SEj] = manual_roc(true,classj);
        
        % Get correlation between ratings
        rN = corr(classi(true==0),classj(true==0),'rows','pairwise');
        rA = corr(classi(true==1),classj(true==1),'rows','pairwise');
        r_avg = (rA+rN)/2;
        
        % get average areas
        A = (Wi+Wj)/2;
            
        
    end
    
end

end


