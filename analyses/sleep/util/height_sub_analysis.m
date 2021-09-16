function levels = height_sub_analysis(pairs_to_plot)

% this function takes an nx2 list of pairs_to_plot, where each row is a
% pair I want to show the stats for, and the members in that row are the
% members of the pair to plot. It returns heights.

% Put in order (1st column before 2nd)
old_pairs_to_plot = pairs_to_plot;
pairs_to_plot(:,1) = min(old_pairs_to_plot,[],2);
pairs_to_plot(:,2) = max(old_pairs_to_plot,[],2);

% Find overlapping groups (groups that would overlap on the x-axis)
levels = nan(size(pairs_to_plot,1),1);
levels(1) = 1;
highest_level = 1;
for i = 2:size(pairs_to_plot)
    
    non_overlap_segment = 0;
    
    % loop over ones already done
    for j = 1:i-1
        % check for overlap
        overlap = find_overlapping_groups(pairs_to_plot([i,j],:));
        
        if overlap == 0
            non_overlap_segment = j;
            break
        end
    end
    
    if non_overlap_segment == 0
        levels(i) = highest_level + 1;
        highest_level = levels(i);
    else
        levels(i) = levels(non_overlap_segment);
    end
    
end


end