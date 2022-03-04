function sz_same_lat =spike_lat_predict_soz_lat(broad_lats,broad_regions,perc_spikes_broad,soz_broad)

sz_same_lat = nan(length(broad_lats),2);

% Loop over spike lateralizations
for i = 1:length(broad_lats)
   
    % Get the indices for the other SPIKE locs (to compare spike rates)
    alt_1 = 3-i;
    
    % get % of spikes in each region
    curr_sp_region = contains(broad_regions,broad_lats{i});
    alt_region_1 = contains(broad_regions,broad_lats{alt_1});
    
    spikes_in_region = nanmean(perc_spikes_broad(curr_sp_region,:),1);
    spikes_outside_region_1 = nanmean(perc_spikes_broad(alt_region_1,:),1);
    

    % Find those patients for whom there are more spikes in this region
    % than the other two
    patients_most_spikes = spikes_in_region > spikes_outside_region_1;
    
    soz_of_these_patients = soz_broad(patients_most_spikes);
    
    % Find ones where SOZ matches current loc
    match = contains(soz_of_these_patients,broad_lats{i});
    
    sz_same_lat(i,:) = [sum(match) sum(patients_most_spikes)];
end


end