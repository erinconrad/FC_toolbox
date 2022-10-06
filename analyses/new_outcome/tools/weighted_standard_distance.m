function SDw = weighted_standard_distance(locs,weights)

n = size(locs,1);

locs(locs>1e6) = nan;

% unweighted version
if isempty(weights)
    weights = ones(n,1);
end

if all(isnan(locs),'all') || all(isnan(weights))
    SDw = nan;
    return
end

%% Calculate weighted mean center
wMean = weighted_mean_center(locs,weights);

nan_rows = any(isnan(locs),2) | (isnan(weights));
L = locs; %L(nan_rows,:) = 0;
W = weights * ones(1,size(L,2)); %W(nan_rows,:) = 0;

L(nan_rows,:) = [];
W(nan_rows,:) = [];
n = length(W);

SDw = sqrt(sum(sum(W.*(L-repmat(wMean,n,1)).^2,1)/sum(W,1),2));

if 0
    hold off
    scatter3(L(:,1),L(:,2),L(:,3),100,W(:,1),'filled');
    hold on
    scatter3(wMean(1),wMean(2),wMean(3),SDw,'k','filled');
    title(sprintf('SD: %1.1f',SDw))
    pause
end

end