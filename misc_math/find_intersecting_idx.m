function [indices,intersecting_set] = find_intersecting_idx(A)

nsets = length(A);
indices = cell(nsets,1);
intersecting_set = A{1};

for i = 2:nsets
    
    intersecting_set = intersect(intersecting_set,A{i});
    
end

for i = 1:nsets
    indices{i} = ismember(A{i},intersecting_set);
end