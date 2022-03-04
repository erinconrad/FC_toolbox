function sz_same_loc = spike_loc_predict_soz_loc(broad_locs,broad_regions,perc_spikes_broad,soz_broad)

sz_same_loc = nan(length(broad_locs),2);
nlocs = length(broad_locs);

% Loop over spike localizations
for i = 1:nlocs
   
    % Get the indices for the other SPIKE locs (to compare spike rates)
    
    alt_1 = mod(i+1,length(broad_locs))+1;
    alt_2 = mod(i,length(broad_locs))+1;
    
    % get % of spikes in each region
    curr_sp_region = contains(broad_regions,broad_locs{i});
    alt_region_1 = contains(broad_regions,broad_locs{alt_1});
    alt_region_2 = contains(broad_regions,broad_locs{alt_2});
    
    spikes_in_region = nanmean(perc_spikes_broad(curr_sp_region,:),1);
    spikes_outside_region_1 = nanmean(perc_spikes_broad(alt_region_1,:),1);
    spikes_outside_region_2 = nanmean(perc_spikes_broad(alt_region_2,:),1);
    

    % Find those patients for whom there are more spikes in this region
    % than the other two
    patients_most_spikes = spikes_in_region > max([spikes_outside_region_1;spikes_outside_region_2],[],1);
    
    soz_of_these_patients = soz_broad(patients_most_spikes);
    
    % Find ones where SOZ matches current loc
    match = contains(soz_of_these_patients,broad_locs{i});
    
    sz_same_loc(i,:) = [sum(match) sum(patients_most_spikes)];
end

end