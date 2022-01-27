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
        
        assert(isequal(lats,repmat({'L';'R'},123,1)))
        
    case 'aal'
        locs{1} = regions{1};
        for i = 2:nregions
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
        
    case 'aal_bernabei'
        
        for i = 1:nregions
            if strcmp(regions{i}(end-1:end),'_R') || strcmp(regions{i}(end-1:end),'_L')
                lats{i} = regions{i}(end);
                locs{i} = regions{i}(1:end-2);
            else
                locs{i} = regions{i};
            end
            
        end
    
end

if 0
    table(regions,locs,lats,loc_nums)
end

end