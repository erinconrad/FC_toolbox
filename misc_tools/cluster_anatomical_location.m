function [loc,lat] = cluster_anatomical_location(anatomy)

loc = cell(length(anatomy),1);
lat = cell(length(anatomy),1);

for ich = 1:length(anatomy)
    curr = anatomy{ich};
    
    if isempty(curr)
        lat{ich} = 'unspecified';
        loc{ich} = 'unspecified';
        continue
    end
    
    
    %% Get laterality
    if contains(curr,'right','IgnoreCase',true) || strcmp(curr(1),'r') || strcmp(curr(1),'R')
        lat{ich} = 'Right';
    end
    
    if contains(curr,'left','IgnoreCase',true) || strcmp(curr(1),'l') || strcmp(curr(1),'L')
        if strcmp(lat{ich},'R')
            error('check laterality');
        else
            lat{ich} = 'Left';
        end
    end
    
    % Make sure I found laterality
    if isempty(lat{ich})
        lat{ich} = 'unspecified'; 
    end
    
    %% Get localization
    if contains(curr,'white','IgnoreCase',true)
        if ~isempty(loc{ich})
            error('check loc');
        end
        loc{ich} = 'white matter';
    end
    
    if contains(curr,'amy','IgnoreCase',true) || ...
            contains(curr,'hipp','IgnoreCase',true) || ...
            contains(curr,'entorhinal','IgnoreCase',true) || ...
            (contains(curr,'temp','IgnoreCase',true) && ...
            (contains(curr,'med','IgnoreCase',true) || ...
            contains(curr,'mes','IgnoreCase',true)))
        if ~isempty(loc{ich})
            error('check loc');
        end
        loc{ich} = 'mesial temporal';
    end

    if (contains(curr,'temp','IgnoreCase',true) && ...
            (~contains(curr,'med','IgnoreCase',true) && ...
            ~contains(curr,'mes','IgnoreCase',true))) || ...
            contains(curr,'planum polare','IgnoreCase',true) || ...
            contains(curr,'fusiform','IgnoreCase',true)
        if ~isempty(loc{ich})
            error('check loc');
        end
        loc{ich} = 'temporal neocortical';
    end
    

    
    % if don't find it, some other checks
    if isempty(loc{ich})
        
        if contains(curr,'gyrus','IgnoreCase',true) || ...
                contains(curr,'cortex','IgnoreCase',true)
            loc{ich} = 'other cortex';
        else
        
            loc{ich} = 'other';
        end
    end
    
    
end

end