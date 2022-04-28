function non_para_roc_comp(models,true)

%{
follows https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-77
compute sd(θ1 - θ2) with N (defaults to 2000) bootstrap replicates. In each replicate r, the original measurements are resampled with replacement; both new ROC curves corresponding to this new sample are built, the resampled AUCs θ1,rand θ2,rand their difference D r = θ1,r- θ2,rare computed. Finally, we compute sd(θ1 - θ2) = sd(D). As Z approximately follows a normal distribution, one or two-tailed p-values are calculated accordingly. This bootstrap test is very flexible and can be applied to AUC, pAUC and smoothed ROC curves.

%}


nboot = 1e3;
nmeas = length(true);



nmodels = length(models);
W = nan(nmodels,1);

for im = 1:nmodels
    classi = models{im};
    
    % Get AUC of model i
    Wi = manual_roc(true,classi);
    W(im) = Wi;
end

D = nan(nmodels,nmodels,nboot);

for ib = 1:nboot
    
    % resample with replacement
    y = randsample(nmeas,nmeas);
    trueb = true(y);
    
    
    
    for im = 1:nmodels
        classi = models{im};
        classi = classi(y);
        Wi = manual_roc(trueb,classi);
              
        for jm = 1:im-1
            classj = models{jm};
            classj = classj(y);
            Wj = manual_roc(trueb,classj);
            
            D(im,jm,ib) = Wi-Wj;
            D(jm,im,ib) = Wi-Wj;
        end
        
    end
end

SD = nan(nmodels,nmodels);
z = nan(nmodels,nmodels);

for im = 1:nmodels
    for jm = 1:im-1
        SD(im,jm) = nanstd(D(im,jm,:));
        SD(jm,im) = nanstd(D(im,jm,:));
        
        z(im,jm) = (Wi-Wj)/SD(im,jm);
        z(jm,im) = (Wi-Wj)/SD(im,jm);
    end
end

end


