function thing_bin = bin_univariate_atlas(thing,atlas_designation,atlas_names)

nregions = length(atlas_names);
thing_bin = nan(nregions,1);

for i = 1:nregions
    % find electrodes in that region
    curr_region = atlas_names{i};
    matching_elecs = strcmp(atlas_designation,curr_region);
    
    % average
    thing_bin(i) = nanmean(thing(matching_elecs));
    
end

end