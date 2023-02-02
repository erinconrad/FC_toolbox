function boxplot_with_points(x,categories,show_stats)

%{
cols = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250];
%}
cols = colormap(gca,lines(length(unique(categories))));

boxplot(x,categories,'colors',cols,'symbol','');
hold on
unique_lats = xticklabels;
nlats = length(unique_lats);
for il = 1:nlats
    curr_lats = strcmp(categories,unique_lats{il}); 
    plot(il + randn(sum(curr_lats),1)*0.05,x(curr_lats),'o','color',cols(il,:))
end
plot(xlim,[0 0],'k--')

if show_stats
    yl = ylim;
    new_y = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];
    ylim(new_y)

    p = kruskalwallis(x,categories,'off');
    bon_p = 0.05/3;
    if p < 0.05
            % do post hoc
            lrp = ranksum(x(strcmp(categories,'left')),x(strcmp(categories,'right')));
            rbp = ranksum(x(strcmp(categories,'right')),x(strcmp(categories,'bilateral')));
            lbp = ranksum(x(strcmp(categories,'left')),x(strcmp(categories,'bilateral')));

            ybar1 = yl(1) + 1.06*(yl(2)-yl(1));
            ytext1 = yl(1) + 1.09*(yl(2)-yl(1));
            ybar2 = yl(1) + 1.18*(yl(2)-yl(1));
            ytext2 = yl(1) + 1.21*(yl(2)-yl(1));
            if lrp < bon_p
                plot([1 2],[ybar1 ybar1],'k-','linewidth',1)
                text(1.5,ytext1,get_asterisks_bonferroni(lrp,3),'horizontalalignment','center','fontsize',15)
            end
            if rbp < bon_p
                plot([2 3],[ybar1 ybar1],'k-','linewidth',1)
                text(2.5,ytext1,get_asterisks_bonferroni(rbp,3),'horizontalalignment','center','fontsize',15)
            end
            if lbp < bon_p
                plot([1 3],[ybar2 ybar2],'k-','linewidth',1)
                text(2,ytext2,get_asterisks_bonferroni(lbp,3),'horizontalalignment','center','fontsize',15)
            end
    else
        ybar1 = yl(1) + 1.06*(yl(2)-yl(1));
        ytext1 = yl(1) + 1.09*(yl(2)-yl(1));
        plot([1 3],[ybar1 ybar1],'k-','linewidth',1)
        text(2,ytext1,'ns','horizontalalignment','center','fontsize',15)

    end

end


end