function elecs_in_parcel = get_elecs_in_parcel(atlas_nums,elecs_atlas)

npts = length(elecs_atlas);
nregions = length(atlas_nums);

elecs_in_parcel = cell(nregions,npts);

for ip = 1:npts
    curr_pt = elecs_atlas{ip};
    for ir = 1:nregions
        elecs_in_parcel{ir,ip} = curr_pt == atlas_nums(ir);
        
    end
end


end