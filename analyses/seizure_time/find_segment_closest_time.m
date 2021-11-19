function find_segment_closest_time(segments,times,szs)

%{
This function takes an tx1 list of times_of_interest and an nxk list of segments and it
outputs an indeterminate length one dimensional vector listing the closest
segment for each time that falls somewhere in the segments.
%}

nsleeps = size(segments,1);
nbins = size(segments,2);
nszs = length(szs);

sz_bins = [];

for is = 1:nsleeps
    curr_sleep = segments(is,:);
    curr_times = times(curr_sleep);
    
    % Loop over seizures
    for iz = 1:nszs
        % see if it falls somewhere in the times
        curr_sz = szs(iz);
        
        if curr_sz >= curr_times(1) && curr_sz <= curr_times(end)
            % find the closest one
            [~,I] = min(abs(curr_times-curr_sz));
            
            sz_bins = [sz_bins;curr_sleep(I)];
            
        end
    end
end

end