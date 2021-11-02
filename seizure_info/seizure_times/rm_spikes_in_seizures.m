function [gdf,n_ictal] = rm_spikes_in_seizures(gdf,fs,run_start,sz_times)

% run start is the time of the run start relative to the file start

% sz_times are the times of the seizures relative to the file start

% gdf(:,2) are the indices into the run of the spikes (need to convert to
% seconds relative to file start)

if isempty(sz_times)
    n_ictal = 0;
    return
end

sp_indices = gdf(:,2);
sp_times = (sp_indices-1)/fs + run_start;

% binary vector indicating, for each spike, if it is in any of the seizures
ictal = (any(sp_times' >= sz_times(:,1) & sp_times' <= sz_times(:,2),1))';

gdf(ictal,:) = [];
n_ictal = sum(ictal);

end