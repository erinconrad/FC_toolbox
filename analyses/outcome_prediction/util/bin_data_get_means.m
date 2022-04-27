function [xbins,ybins] = bin_data_get_means(x,y,percentile)



%% Define bins
bins = min(x):(percentile*(max(x)-min(x))):max(x);
nbins = length(bins);
xbins = nan(nbins-1,1);
ybins = nan(nbins-1,1);
for ib = 1:nbins-1
    xbins(ib) = mean([bins(ib) bins(ib+1)]);
    ybins(ib) = nanmean(y(x>=bins(ib) & x < bins(ib+1)));
end


end