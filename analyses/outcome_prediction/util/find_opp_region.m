function opp_loc = find_opp_region(ir,curr_loc_i,locs,lats)

 % find the other one with same name
all_locs_same = find(strcmp(locs,curr_loc_i));
opp_loc = all_locs_same(~ismember(all_locs_same,ir));

if isempty(opp_loc)
    assert(isempty(lats{ir}))
    return; 
end

% confirm it's opposite side
assert(~strcmp(lats{ir},lats{opp_loc}))

end