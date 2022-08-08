function comb = homogenize_soz_locs_lats(loc,lat,which_localization)

switch which_localization
    case 'fine'
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
        
    case 'broad'
        if contains(loc,'temporal')
            loc = 'temporal';
        elseif strcmp(loc,'other cortex')
            loc = 'other cortex';
        elseif strcmp(loc,'multifocal') || contains(loc,'diffuse')
            loc = 'multifocal/diffuse';
        else
            loc = [];
        end
        
    case 'broad_soz'
        if contains(loc,'temporal')
            loc = 'temporal';
        elseif strcmp(loc,'other cortex')
            loc = 'other cortex';
        elseif strcmp(loc,'multifocal') || contains(loc,'diffuse')
            loc = 'multifocal/diffuse';
        else
            loc = [];
        end

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

if strcmp(which_localization,'broad') || strcmp(which_localization,'broad_soz')
    if contains(comb,'diffuse')
        comb = 'bilateral/diffuse';
    end
end

end