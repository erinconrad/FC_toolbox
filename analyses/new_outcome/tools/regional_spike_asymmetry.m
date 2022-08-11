function AI = regional_spike_asymmetry(spike_rates)
%{
This function takes a distribution of spike rates across regions and
calculates a measure of asymmetry (how different the distribution is from a
uniform distribution).

Inputs:
- spike_rates: an Nx1 vector of average spike rates in each region, where N is
number of regions. 

Think about how to handle nans. My thought is that they should not be
counted in the average or in the difference part. Doing this means I should
NOT have to control for the number of electrodes in each region.

This assumes under the null a uniform distribution of spikes, which of
course may not be true (there may be regional variation in spikes that
occurs by default).
%}

nregions = length(spike_rates);

numerator = 0;
denominator = 0;

% If only 1 region without nans, make whole thing nans
if sum(~isnan(spike_rates)) <=1
    AI = nan; 
    return;
end

for i = 1:nregions
    if isnan(spike_rates(i)), continue; end
        
    % Add the absolute value of the difference between the spike rates in
    % this region and the mean across the non-nan spikes
    numerator = numerator + abs(spike_rates(i)-mean(spike_rates(~isnan(spike_rates))));
    
    % add the spike rates in this region
    denominator = denominator + spike_rates(i);
end

AI = numerator/denominator;

end