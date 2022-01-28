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
        
    case 'brainnetome'
        
        for i = 1:length(names)
            curr = names{i};
            
            if contains(curr,'Tha') || contains(curr,'BG')
            elseif contains(curr,'Amyg') || contains(curr,'Hipp') || ...
                    contains(curr,'PhG')
                if strcmp(curr(end-1:end),'_L')
                    locs{i} = 'left mesial temporal';
                elseif strcmp(curr(end-1:end),'_R')
                    locs{i} = 'right mesial temporal';
                end
            elseif contains(curr,'STG') || contains(curr,'MTG') || ...
                    contains(curr,'ITG') || contains(curr,'FuG') || ...
                    contains(curr,'pSTS')
                if strcmp(curr(end-1:end),'_L')
                    locs{i} = 'left temporal neocortical';
                elseif strcmp(curr(end-1:end),'_R')
                    locs{i} = 'right temporal neocortical';
                end
            else
                if strcmp(curr(end-1:end),'_L')
                    locs{i} = 'left other cortex';
                elseif strcmp(curr(end-1:end),'_R')
                    locs{i} = 'right other cortex';
                end
            end
        
        end
                 
end

if 0
    table(names,locs)
end

end