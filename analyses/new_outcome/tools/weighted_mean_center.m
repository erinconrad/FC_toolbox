function weightedmean = weighted_mean_center(L,W)

nan_rows = any(isnan(L),2) | (isnan(W));
Lf = L; Lf(nan_rows,:) = 0;
Wf= W * ones(1,size(L,2)); Wf(nan_rows,:) = 0;
weightedmean = sum(Lf .* Wf, 1) ./ sum(Wf, 1);

if 0
    figure
    scatter3(L(:,1),L(:,2),L(:,3),100,W,'filled');
    hold on
    scatter3(weightedmean(1),weightedmean(2),weightedmean(3),200,'k','filled');
    
end

end