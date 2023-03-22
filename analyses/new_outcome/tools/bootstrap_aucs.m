function aucs = bootstrap_aucs(class,scores,pos_class,nb)

aucs = nan(nb,1);
npts = length(scores);
for ib = 1:nb
    % sample the patients with replacement. Put in a while loop because I
    % want to resample if I only get one class
    while 1
        y = randsample(npts,npts,true); % true means sample with replacement
        temp_class = class(y);
        temp_scores = scores(y);
        if length(unique(temp_class)) == 1
            continue % try again
        else
            [~,~,~,auc] = perfcurve(temp_class,temp_scores,pos_class);
            aucs(ib) = auc;
            break
        end
    end

end


end