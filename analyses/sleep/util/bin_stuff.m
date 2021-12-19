function [bin_counts,closest_bins] = bin_stuff(stuff,bins)

bin_counts = zeros(size(bins));
closest_bins = nan(size(stuff));
for i = 1:length(stuff)
    
    % find the closest without going over
    time_diff = stuff(i) - bins;
    time_diff(time_diff >0) = -inf;
    [~,closest] = max(time_diff);
    
    if sum(time_diff == -inf) == length(time_diff)
        closest_bins(i) = nan;
    else
    
        closest_bins(i) = closest;
        bin_counts(closest) = bin_counts(closest) + 1;
    
    end
        
end

end