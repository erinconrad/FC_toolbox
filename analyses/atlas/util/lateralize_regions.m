function [locs,lats,loc_nums] = lateralize_regions(regions,atlas)

nregions = length(regions);
locs = cell(nregions,1);
lats = cell(nregions,1);
loc_nums = nan(nregions,1);

loc_nums(1) = 1;

switch atlas
    
    case 'brainnetome'
    
        for i = 1:nregions
            lats{i} = regions{i}(end);
            locs{i} = regions{i}(1:end-2);
            if i > 1
                % if it's the same loc as the last, make number same
                if isequal(locs{i},locs{i-1})
                    loc_nums(i) = loc_nums(i-1);
                else
                    % advance it by one
                    loc_nums(i) = loc_nums(i-1) + 1;
                end  
            end
            
        end
    
end


end