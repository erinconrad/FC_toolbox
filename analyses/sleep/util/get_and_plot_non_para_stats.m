function [p,post_hoc_p,which_groups]=get_and_plot_non_para_stats(yl,data)

[p,post_hoc_p,which_groups] = non_parametric_stats(data);
if p > 0.05
    pairs_to_plot = [];
else
    pairs_to_plot = which_groups(post_hoc_p < 0.05/size(which_groups,1),:);
    post_hoc_p_to_plot = post_hoc_p(post_hoc_p < 0.05/size(which_groups,1));
end

heights = get_heights(yl,pairs_to_plot);

if p > 0.05
    plot([1 size(data,1)],...
        [heights(size(heights,1)-1,1) heights(size(heights,1)-1,1)],'k',...
        'linewidth',2);
    text(mean([1 size(data,1)]),heights(size(heights,1)-1,2),...
        'ns','fontsize',15,'horizontalalignment','center')
else
    for k = 1:size(pairs_to_plot,1)
        plot([pairs_to_plot(k,1)+0.1 pairs_to_plot(k,2)-0.1],[heights(k,1) heights(k,1)],'k-',...
            'linewidth',2)
        hold on
        text(mean(pairs_to_plot(k,:)),heights(k,2),...
            get_asterisks(post_hoc_p_to_plot(k),size(which_groups,1)),...
            'fontsize',15,'horizontalalignment','center')
    end
end

ylim([yl(1) heights(end,2)])

end