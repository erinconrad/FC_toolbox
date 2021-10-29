function plot_and_stats_change(thing,xtext,ytext,p_or_unp)

% thing is ngroups x npts x 2 (where 2 is the thing we're measuring change
% over)
ngroups = size(thing,1);
npts = size(thing,2);

do_rel = 0;

change = nan(ngroups,npts);

for ig = 1:ngroups
    curr = squeeze(thing(ig,:,:));
    if do_rel
        change(ig,:) = (curr(:,2) - curr(:,1))./abs(curr(:,1));
    else
        change(ig,:) = curr(:,2) - curr(:,1);
    end   
end

% plot and stats
plot_paired_data(change,xtext,ytext,'paired')


end