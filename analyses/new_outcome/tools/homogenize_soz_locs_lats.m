function comb = homogenize_soz_locs_lats(loc,lat)


if strcmp(loc,'mesial temporal')
    loc = 'mesial temporal';
elseif strcmp(loc,'other cortex')
    loc = 'other cortex';
elseif strcmp(loc,'temporal neocortical')
    loc = 'temporal neocortical';
elseif strcmp(loc,'multifocal') || contains(loc,'diffuse')
    loc = 'multifocal/diffuse';
else
    loc = [];
end

if strcmp(lat,'left')
    lat = 'left';
elseif strcmp(lat,'right')
    lat = 'right';
elseif strcmp(lat,'bilateral') || contains(lat,'diffuse')
    lat = 'bilateral/diffuse';
else
    lat = [];
end
       

%% Combine them
if isempty(loc) && isempty(lat)
    comb = ' ';
elseif strcmp(lat,'bilateral/diffuse') || (isempty(lat) && strcmp(loc,'multifocal/diffuse'))
    
    comb = 'bilateral/diffuse';
elseif isempty(loc)
    comb = sprintf('%s multifocal/diffuse',lat);
else
    comb = sprintf('%s %s',lat,loc);
end


end