function [loc,lat] = seizure_localization_parser(soz_loc,soz_lat)

if isempty(soz_loc)
    loc = nan;
elseif contains(soz_loc,'temporal')
    loc = 'temporal';
else
    loc = 'other';
end

lat = soz_lat;

end