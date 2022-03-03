function out_locs = bin_into_temp_extra(locs)

out_locs = cell(size(locs));

for i = 1:length(locs)
    curr = locs{i};
    
    if isempty(curr), continue; end
    
    if contains(curr,'temporal')
        out_locs{i} = 'temporal';
    elseif contains(curr,'other cortex')
        out_locs{i} = 'other cortex';
    end
end

end