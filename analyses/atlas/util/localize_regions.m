function locs = localize_regions(names,atlas)

locs = cell(length(names),1);

switch atlas
    
    case 'aal_bernabei'
        for i = 1:length(names)
            curr = names{i};
            
            if contains(curr,'Hippocamp') || contains(curr,'Amygdala')
                if strcmp(curr(end-1:end),'_L')
                    locs{i} = 'left mesial temporal';
                elseif strcmp(curr(end-1:end),'_R')
                    locs{i} = 'right mesial temporal';
                end
            elseif contains(curr,'Fusiform') || contains(curr,'Heschl') || ...
                    contains(curr,'Temporal')
                if strcmp(curr(end-1:end),'_L')
                    locs{i} = 'left temporal neocortical';
                elseif strcmp(curr(end-1:end),'_R')
                    locs{i} = 'right temporal neocortical';
                end
                
            elseif contains(curr,'White_Matter') || contains(curr,'Vermis') || ...
                    contains(curr,'Cerebelum') || contains(curr,'Thalamus') || ...
                    contains(curr,'Putamen') || contains(curr,'Caudate') || ...
                    contains(curr,'Pallidum')
                %locs{i} = 'other';
            else
                if strcmp(curr(end-1:end),'_L')
                    locs{i} = 'left other cortex';
                elseif strcmp(curr(end-1:end),'_R')
                    locs{i} = 'right other cortex';
                end
                
            end
                
        end
        
        if 1
            table(names,locs)
        end
        
end

end