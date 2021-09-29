function indices = find_pts_with_scalp(pt)

names = {};
indices = [];

for p = 1:length(pt)
    
    if isempty(pt(p).ieeg.file(1).chLabels)
        continue
    end
    elecs = pt(p).ieeg.file(1).chLabels;
    if any(ismember(elecs,'CZ')) && any(ismember(elecs,'C3')) && any(ismember(elecs,'C4'))
        names = [names;pt(p).name];
        indices = [indices;p];
    end
    
end


end