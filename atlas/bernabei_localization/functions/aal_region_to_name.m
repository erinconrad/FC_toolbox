function [names,regions] = aal_region_to_name(filepath,regions)

T = readtable(filepath,'ReadVariableNames',false);
names = cell(length(regions),1);

if isempty(regions)
    regions = T.Var3;
end

for i = 1:length(regions)
    row = T.Var3 == regions(i);
    assert(sum(row) <= 1)
    if sum(row) == 1
        names{i} = T.Var2{row};
    end
    
end

end