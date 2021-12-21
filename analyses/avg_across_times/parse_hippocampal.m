function out = parse_hippocampal(out)

for i = 1:length(out)
    if ~isempty(out(i).anatomy)
        [left_hipp,right_hipp] = add_hippocampal(out(i).anatomy);
        out(i).left_right_hipp = [left_hipp right_hipp];
    end
    
end


end




function [left_hipp,right_hipp] = add_hippocampal(labels)

right_hipp = zeros(length(labels),1);
left_hipp = zeros(length(labels),1);

for i = 1:length(labels)
    if isempty(labels{i}), continue; end
    if contains(labels{i},'hippocampus','IgnoreCase',true)
        if contains(labels{i},'Left','IgnoreCase',true)
            left_hipp(i) = 1;
        elseif contains(labels{i},'right','IgnoreCase',true)
            right_hipp(i) = 1;
        end
    end
    
end



end

