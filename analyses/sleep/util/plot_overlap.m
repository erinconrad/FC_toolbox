function plot_overlap(overlap,xtext,ytext,nb)

ci = bootci(nb,@(x) mean(x),overlap);
errorbar(1:size(overlap,2),mean(overlap,1),...
    mean(overlap,1)-ci(1,:),...
    ci(2,:)-mean(overlap,1),'o','linewidth',2)
xlim([0 3])
xticks([1 2])
xticklabels(xtext)
ylabel(ytext)

end