function interaction_plot_and_stats(thing,xtext,ytext)

thing = permute(thing,[3 2 1]);

% loop groups
for i = 1:size(thing,1)
    errorbar([i-0.2 i+0.2],squeeze(nanmean(thing(i,:,:),2)),...
        squeeze(nanstd(thing(i,:,:),[],2)),'o','linewidth',2)
    hold on
end
xticks(1:size(thing,1))
xticklabels(xtext)
ylabel(ytext)

end