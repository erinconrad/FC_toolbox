function bins = bin_seizure_times(sz_times,times,surround,rm_cluster,surround_secs)

%{
Take a list of times and and a list of seizure times and output a list of
bins, one bin for each seizure, with surround x 2 elements in each bin
listing the indices surrounding the seizure. I can either include or
exclude clustered seizures
%}

%% If I choose, remove clustered seizures
nszs = size(sz_times,1);
if rm_cluster
    sz_to_rm = zeros(nszs,1);
    all_sz_idx = (1:nszs)';
    for s = 1:nszs
        loo = all_sz_idx ~= s;
        if any(abs(sz_times(s,1) - sz_times(loo,1)) < surround_secs)
            sz_to_rm(s) = 1;
        end
    end
    sz_times(logical(sz_to_rm),:) = [];
    nszs = size(sz_times,1);
end


%% Loop over remaining szs
bins = nan(nszs,surround*2);
for s = 1:nszs
    % find the closest time
    [~,closest] = min(abs(times-sz_times(s,1)));
    
    bins(s,:) = closest - surround : closest + surround-1;
end

end