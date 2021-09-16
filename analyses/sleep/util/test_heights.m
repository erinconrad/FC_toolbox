function test_heights

things = [1 2 3 4 5 6];
sig_comp = [1 2;2 3;1 3];

plot(things,ones(length(things),1),'o','markersize',20,'linewidth',2)
hold on
yl = ylim;
heights = get_heights(yl,sig_comp);

% main
plot([things(1) things(end)],...
    [heights(size(heights,1)-1,1) heights(size(heights,1)-1,1)],'k-')
hold on
text(mean(things),heights(size(heights,1)-1,2),get_asterisks(0.01,1),'fontsize',20)

for i = 1:size(sig_comp,1)
    plot([sig_comp(i,1)+0.1 sig_comp(i,2)-0.1],[heights(i,1) heights(i,1)],'k-')
    hold on
    text(mean(sig_comp(i,:)),heights(i,2),get_asterisks(0.01,1),'fontsize',20)
end

ylim([yl(1) heights(end,2)])

end