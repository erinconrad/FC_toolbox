function data = remove_if_nans_surround(data,span_to_look,max_nans)

nruns = size(data,2);

% loop over runs
for r = span_to_look+1:nruns-span_to_look
    
    % Take span
    span = data(:,r-span_to_look:r+span_to_look);
    
    % Find the rows for which there are more than max_nans nans
    too_many_nans = sum(isnan(span),2) > max_nans;
    
    % turn the middle to a nan as well (note this may have a recursive
    % effect
    data(too_many_nans,r) = nan;
    
end
    

end