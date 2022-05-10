function [new_loc,new_lat] = bin_manual_locs(loc,lat)

if isempty(loc) || isempty(lat)
    new_loc = 'no seizures';
    new_lat = 'no seizures';
    return
end

switch loc
    case 'mesial temporal'
        new_loc = 'mesial temporal';
    case 'temporal neocortical'
        new_loc = 'temporal neocortical';
    case 'other cortex'
        new_loc = 'other cortex';
    otherwise
        new_loc = 'other';
end

switch lat
    case 'left'
        new_lat = 'left';
    case 'right'
        new_lat = 'right';
    case 'bilateral'
        new_lat = 'bilateral';
    case 'diffuse'
        new_lat = 'bilateral';
    otherwise
        error('what');
end
    

end