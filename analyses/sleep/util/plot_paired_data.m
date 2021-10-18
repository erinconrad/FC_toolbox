function plot_paired_data(data,xlabels,ytext)

%errorbar(nanmean(data,2),nanstd(data,[],2),'o','linewidth',2);
pr = prctile(data,[25,75],2);
errorbar(1:size(data,1),nanmedian(data,2),...
    nanmedian(data,2)-pr(:,1),pr(:,2)-nanmedian(data,2),'o','linewidth',2);
hold on
xlim([0 size(data,1)+1])
xticks(1:size(data,1))
xticklabels(xlabels)
ylabel(ytext)
yl = ylim;
[p,post_hoc_p,which_groups]=get_and_plot_non_para_stats(yl,data);
set(gca,'fontsize',15)

end