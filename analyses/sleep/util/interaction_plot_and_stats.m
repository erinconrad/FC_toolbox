function interaction_plot_and_stats(thing,legend_text,ytext,xtext,show_leg)

%thing = permute(thing,[3 2 1]);

% loop groups
lp = nan(size(thing,1),1);
xticklocs = [];
xtextreap = {};
for i = 1:size(thing,1)
    %plot([i-0.2 i+0.2],
    %
    lp(i) = errorbar([i-0.2 i+0.2],squeeze(nanmedian(thing(i,:,:),2)),...
        squeeze(iqr(thing(i,:,:),2)),'o','linewidth',2);
    %{
    lp(i) = errorbar([i-0.2 i+0.2],squeeze(nanmean(thing(i,:,:),2)),...
        squeeze(nanstd(thing(i,:,:),[],2)),'o','linewidth',2);
    %}
    xticklocs = [xticklocs,i-0.2,i+0.2];
    xtextreap = [xtextreap,xtext{1},xtext{2}];
    %}
    hold on
end
xticks(xticklocs)
xticklabels(xtextreap)
xlim([0 size(thing,1)+1])
xtickangle(45)
ylabel(ytext)
if show_leg
    legend(lp,legend_text,'location','northeastoutside','fontsize',15);
end
set(gca,'fontsize',15)

end