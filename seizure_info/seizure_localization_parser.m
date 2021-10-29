function [loc,lat] = seizure_localization_parser(soz_loc,soz_lat)

if contains(soz_loc,'temporal')
    loc = 'temporal';
else
    loc = 'other';
end

lat = soz_lat;

end