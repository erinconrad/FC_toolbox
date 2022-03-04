function same_loc = soz_loc_predict_spike_loc(soz_broad,broad_locs,perc_spikes_broad,broad_regions)

same_loc = nan(length(broad_locs),2);

for i = 1:length(broad_locs)

    % Get patients with that sz loc
    pts_with_sz_loc = contains(soz_broad,broad_locs{i});
    
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
    
    any_isnot_nan = ~any(~isnan([spikes_in_region;spikes_outside_region_1;spikes_outside_region_2]),1);
    
    % For which patients is the % of spikes highest in this region?
    most_spikes_in_region = (spikes_in_region > max([spikes_outside_region_1;spikes_outside_region_2],[],1));
    
    % % of patient with this sz loc meeting this criteria
    same_loc(i,:) = [sum(most_spikes_in_region(pts_with_sz_loc & ~any_isnot_nan')),sum(pts_with_sz_loc & ~any_isnot_nan')];
    
    if 0
        table(soz_broad(pts_with_sz_loc& ~any_isnot_nan'),...
            (spikes_in_region(pts_with_sz_loc& ~any_isnot_nan'))',...
            (spikes_outside_region_1(pts_with_sz_loc& ~any_isnot_nan'))',...
            (spikes_outside_region_2(pts_with_sz_loc& ~any_isnot_nan'))')
        
        table(soz_broad(pts_with_sz_loc),...
            (spikes_in_region(pts_with_sz_loc))',...
            (spikes_outside_region_1(pts_with_sz_loc))',...
            (spikes_outside_region_2(pts_with_sz_loc))')
        
    end
    
end


end