function [pval,all_overlap,real_overlap] = bootstrap_lr_features(T,features,nb)

nfeatures = length(features);
all_overlap = nan(nb,nfeatures);
response = T.soz_lats;

% do it once for real
feature_p_val = nan(nfeatures,2);
feature_eta2= nan(nfeatures,2);
for i = 1:nfeatures
    curr_feat = T.(features{i});

    dT = cohenD(curr_feat(strcmp(response,'left')),curr_feat(strcmp(response,'bilateral')));
    feature_eta2(i,1) = dT;
    dT = cohenD(curr_feat(strcmp(response,'right')),curr_feat(strcmp(response,'bilateral')));
    feature_eta2(i,2) = dT;

    feature_p_val(i,1) = ranksum(curr_feat(strcmp(response,'left')),curr_feat(strcmp(response,'bilateral')));
    feature_p_val(i,2) = ranksum(curr_feat(strcmp(response,'right')),curr_feat(strcmp(response,'bilateral')));
end


[~,I1] = sort(abs(feature_eta2(:,1)),'descend'); [~,I2] = sort(abs(feature_eta2(:,2)),'descend'); 

real_overlap = nan(nfeatures,1);
for k = 1:nfeatures
    real_overlap(k) = length(intersect(I1(1:k),I2(1:k)));
end



% loop over bootstrap iterations
for ib = 1:nb
    if mod(ib,100) == 0
        fprintf('\nDoing %d of %d\n',ib,nb)
    end
    % Shuffle left and right labels
    left = find(strcmp(response,'left')); nleft = length(left);
    right = find(strcmp(response,'right')); nright = length(right);
    bilateral = find(strcmp(response,'bilateral'));
    left_or_right = [left;right]; % the first nleft of these are left and the rest are right

    % randomly permute these
    p = randperm(nleft+nright);
    shuffled = left_or_right(p);
    fake_left = shuffled(1:nleft);
    fake_right = shuffled(nleft+1:end); assert(length(fake_right) == nright)
    
    % get the effect sizes
    feature_eta2= nan(nfeatures,2);
    for i = 1:nfeatures
        curr_feat = T.(features{i});
        
        
        dT = cohenD(curr_feat(fake_left),curr_feat(bilateral));
        feature_eta2(i,1) = dT;
        dT = cohenD(curr_feat(fake_right),curr_feat(bilateral));
        feature_eta2(i,2) = dT;
    end

    
    % get top k effect sizes
    [~,I1] = sort(abs(feature_eta2(:,1)),'descend'); [~,I2] = sort(abs(feature_eta2(:,2)),'descend'); 
    
    % calculate overlap
    for k = 1:nfeatures
        overlap = length(intersect(I1(1:k),I2(1:k)));
        all_overlap(ib,k) = overlap;
    end
    

    
end

% can I make this two tailed instead?
pval = nan(nfeatures,1);
for k = 1:nfeatures
    pval(k) = (sum(all_overlap(:,k) <= real_overlap(k)) + 1)/(nb+1);
end

if 0
figure
plot(sort(all_overlap),'o')
hold on
plot(xlim,[real_overlap real_overlap],'-')
one_tailed_p = (sum(all_overlap<=real_overlap)+1)/(nb+1);
end


end

