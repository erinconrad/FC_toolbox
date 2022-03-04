function same_lat = soz_lat_predict_spike_lat(soz_broad,broad_lats,perc_spikes_broad,broad_regions)

same_lat = nan(length(broad_lats),2);
for i = 1:length(broad_lats)

    % Get patients with that sz loc
    pts_with_sz_loc = contains(soz_broad,broad_lats{i});
    
    % Get the indices for the other SPIKE locs (to compare spike rates)
    alt_1 = 3-i;
    
    % get % of spikes in each region
    curr_sp_region = contains(broad_regions,broad_lats{i});
    alt_region_1 = contains(broad_regions,broad_lats{alt_1});
    
    spikes_in_region = nanmean(perc_spikes_broad(curr_sp_region,:),1);
    spikes_outside_region_1 = nanmean(perc_spikes_broad(alt_region_1,:),1);
    
    % For which patients is the % of spikes highest in this region?
    most_spikes_in_region = (spikes_in_region > spikes_outside_region_1);
    
    any_nan = any(isnan([spikes_in_region;spikes_outside_region_1]),1);
    
    
    if 0
        table(soz_broad(pts_with_sz_loc& ~any_nan'),(spikes_in_region(pts_with_sz_loc& ~any_nan'))',(spikes_outside_region_1(pts_with_sz_loc& ~any_nan'))')
        
    end
    
    % % of patient with this sz loc meeting this criteria
    same_lat(i,:) = [sum(most_spikes_in_region(pts_with_sz_loc & ~any_nan')),sum(pts_with_sz_loc& ~any_nan')];
    
end


end