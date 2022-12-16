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

   if 0
        la1 = find(strcmp(pt(i).elecs(1).elec_names,'LA1'));
        ra1 = find(strcmp(pt(i).elecs(1).elec_names,'RA1'));
        la1_ra1_dist = vecnorm(locs(la1)-locs(ra1));

        lb1 = find(strcmp(pt(i).elecs(1).elec_names,'LB1'));
        rb1 = find(strcmp(pt(i).elecs(1).elec_names,'RB1'));
        lb1_rb1_dist = vecnorm(locs(lb1)-locs(rb1));

        la1_lb1_dist = vecnorm(locs(lb1)-locs(la1));
        lb2 = find(strcmp(pt(i).elecs(1).elec_names,'LB2'));
        lb1_lb2_dist = vecnorm(locs(lb1)-locs(lb2));
   end

end



end