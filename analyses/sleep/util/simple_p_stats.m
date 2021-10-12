function simple_p_stats(yl,p,xpoints)

ybar = yl(1) + 1.05*(yl(2)-yl(1));
ytext = yl(1) + 1.1*(yl(2) - yl(1));
ylnew = [yl(1) yl(1) + 1.15*(yl(2)-yl(1))];

plot([xpoints(1) xpoints(2)],[ybar ybar],'k-','linewidth',2);
text(mean([xpoints(1) xpoints(2)]),ytext,get_asterisks(p,1),...
    'horizontalalignment','center','fontsize',15);
ylim(ylnew)

end