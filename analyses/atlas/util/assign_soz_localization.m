function [out,full_broad] = assign_soz_localization(sozs,broad_regions)

full_broad = [broad_regions;'left broad';'right broad';'bilateral'];
out = cell(size(sozs));

for i = 1:length(sozs)
    curr = sozs{i};
    
    match = strcmp(curr,full_broad);
    
    % if no match
    if sum(match) == 0
        if contains(curr,'left')
            match(7) = 1;
        elseif contains(curr,'right')
            match(8) = 1;
        elseif contains(curr,'bilateral') || strcmp(curr,'diffuse diffuse')
            match(9) = 1;
        end
    end
    
    if sum(match) == 0 % should only happen if no epilepsy localization
        assert(isempty(curr)||strcmp(curr,' '))
        continue;
    end
    
    assert(sum(match)==1)
    
    out{i} = full_broad{match};
    
end


end