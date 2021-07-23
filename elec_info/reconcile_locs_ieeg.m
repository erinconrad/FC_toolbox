function [ieeg_locs,ieeg_anatomy] = reconcile_locs_ieeg(ieeg_names,loc_names,locs,anatomy)

% initialize ieeg_locs and anatomy
nieeg = length(ieeg_names);
ieeg_locs = nan(nieeg,3);
ieeg_anatomy = cell(nieeg,1);

% Loop through ieeg names
for i = 1:nieeg
    
    % curr_name
    curr_name = ieeg_names{i};
    
    % find matching loc name
    matching_loc_name = strcmp(curr_name,loc_names);
    
    if sum(matching_loc_name) == 0
        % no match, don't give it a location
        continue;
    end
    
    if sum(matching_loc_name) > 1
        error('why are there multiple matching names?');
    end
    
    % get the matching loc
    matching_loc = locs(matching_loc_name,:);
    
    % matching anatomy
    if strcmp(class(anatomy),'double')
        matching_anatomy = [];
    else
        matching_anatomy = anatomy{matching_loc_name};
    end
    
    % fill up output array and cell
    ieeg_locs(i,:) = matching_loc;
    ieeg_anatomy{i} = matching_anatomy;
    
end

end