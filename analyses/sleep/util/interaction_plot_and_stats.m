function interaction_plot_and_stats(thing,legend_text,ytext,xtext,show_leg)

%thing = permute(thing,[3 2 1]);

% loop groups
lp = nan(size(thing,1),1);
xticklocs = [];
xtextreap = {};
for i = 1:size(thing,1)
    %plot([i-0.2 i+0.2],
    %
    pr = squeeze(prctile(thing(i,:,:),[25,75],2));
    lp(i) = errorbar([i-0.2 i+0.2],squeeze(nanmedian(thing(i,:,:),2)),...
        squeeze(nanmedian(thing(i,:,:),2))-pr(:,1),...
        pr(:,2)-squeeze(nanmedian(thing(i,:,:),2)),'o','linewidth',2);
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
set(gca,'fontsize',15)

all_ps = stratification_analysis(thing);

yl = ylim;
ybar = yl(1) + 1.1*(yl(2)-yl(1));
ytext = yl(1) + 1.2*(yl(2)-yl(1));
newyl = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];

for i = 1:size(thing,1)
    plot([i-0.2 i+0.2],[ybar ybar],'k-','linewidth',2)
    text(i,ytext,get_asterisks(all_ps(i),length(all_ps)),...
        'horizontalalignment','center','fontsize',15);
end

ylim(newyl)
if show_leg
    legend(lp,legend_text,'location','northeastoutside','fontsize',15);
end
end