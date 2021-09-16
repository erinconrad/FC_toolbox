function p = find_matching_pt(name,pt)

p = nan;
for i = 1:length(pt)
    if strcmp(name,pt(i).name)
        p = i;
        break
    end
end

end