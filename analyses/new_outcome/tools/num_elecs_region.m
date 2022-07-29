function num_elecs = num_elecs_region(unique_regions,elec_broad_names)

nregions = length(unique_regions);
num_elecs = nan(nregions,1);
for i = 1:nregions
    curr_region = unique_regions{i};
    num_elecs(i) = sum(strcmp(elec_broad_names,curr_region));
    
end

end