function plot_elec_locs(pt,start)

npts = length(pt);
figure
for i = start%:npts
    if ~isfield(pt(i),'elecs'), continue; end
    if ~isfield(pt(i).elecs(1),'locs'), continue; end
    if isempty(pt(i).elecs(1).locs), continue; end
    locs = pt(i).elecs_native(1).locs;
    locs(any(abs(locs)>1e10,2),:) = [];
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'filled')
    text(locs(:,1),locs(:,2),locs(:,3),pt(i).elecs(1).elec_names)
    title(pt(i).name)
   % pause

end



end