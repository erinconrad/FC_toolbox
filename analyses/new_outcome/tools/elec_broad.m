function elec_broad_names = elec_broad(atlas,atlas_names,broad_names)

nelecs = length(atlas);
elec_broad_names = cell(nelecs,1);
for i = 1:nelecs
    curr_atlas_name = atlas{i};
    atlas_index = strcmp(curr_atlas_name,atlas_names);
    if sum(atlas_index) == 0
        elec_broad_names{i} = [];
    else
        elec_broad_names{i} = broad_names{atlas_index};
    end
    
end

end