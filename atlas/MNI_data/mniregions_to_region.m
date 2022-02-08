function regions = mniregions_to_region(region_name)

regions = cell(length(region_name),1);
for i = 1:length(region_name)
    curr = region_name{i};
    if strcmp(curr,'Amygdala') || strcmp(curr,'Hippocampus')
        regions{i} = 'mesial temporal';
    elseif contains(curr,'Fusiform') || contains(curr,'emporal')
        regions{i} = 'temporal neocortical';
    elseif strcmp(curr,'White matter')
    else
        regions{i} = 'other cortex';
    end
    
    
end

if 0
    table(region_name,regions)
end

end