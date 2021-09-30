function indices = find_pts_with_scalp(pt)

names = {};
indices = [];

for p = 1:length(pt)
    
    if isempty(pt(p).ieeg.file(1).chLabels)
        continue
    end
    elecs = pt(p).ieeg.file(1).chLabels;
    if (any(ismember(elecs,'C3')) && any(ismember(elecs,'F3')) && any(ismember(elecs,'T3'))) || ...
            (any(ismember(elecs,'C4')) && any(ismember(elecs,'F4')) && any(ismember(elecs,'T4')))
        names = [names;pt(p).name];
        indices = [indices;p];
    end
    
end


end