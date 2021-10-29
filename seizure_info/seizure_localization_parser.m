function [loc,lat] = seizure_localization_parser(soz_loc,soz_lat)

if strcmp(soz_loc,'mesial temporal')
    loc = 'mesial temporal';
else
    loc = 'other';
end

lat = soz_lat;

end