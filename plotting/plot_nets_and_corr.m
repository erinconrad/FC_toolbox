function plot_nets_and_corr(net1,net2)

% Expand if flat
[m,n] = size(net1);
smallest_dim = min([m,n]);
if smallest_dim == 1
    net1 = wrap_or_unwrap_adjacency_fc_toolbox(net1);
    net2 = wrap_or_unwrap_adjacency_fc_toolbox(net2);
end
    

figure
tiledlayout(2,2)
nexttile
turn_nans_gray(net1)

nexttile
turn_nans_gray(net2)

nexttile([1 2])
ns1 = nansum(net1,2);
ns2 = nansum(net2,2);
[r,pval] = corr(ns1,ns2,'rows','pairwise');
plot(ns1,ns2,'o','linewidth',2);
yl = ylim;
xl = xlim;
text(xl(2),yl(2),sprintf('r = %1.2f, p = %1.3f',r,pval),...
    'horizontalalignment','right','verticalalignment','top');


end