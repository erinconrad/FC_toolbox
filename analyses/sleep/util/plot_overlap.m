function plot_overlap(overlap,xtext,ytext,nb,ps)

%% Get stats
pval = nan(size(overlap,2),1);
for i = 1:size(overlap,2)
    [pval(i),exp] = binomial_generalized(ps,sum(overlap(:,i)));
end
pval_diff = binomial_generalized_difference(ps,sum(diff(overlap,1,2)));

ci = bootci(nb,@(x) mean(x),overlap);
%{
errorbar(1:size(overlap,2),mean(overlap,1),...
    mean(overlap,1)-ci(1,:),...
    ci(2,:)-mean(overlap,1),'o','linewidth',2)
%}
plot(mean(overlap,1),'o','linewidth',2,'markersize',10)
hold on
plot(xlim,[exp exp],'k--','linewidth',2)
xlim([0 3])

yl = ylim;
yast1 = yl(1) + 1.05*(yl(2)-yl(1));
ybar = yl(1) + 1.15*(yl(2)-yl(1));
yast2 = yl(1) + 1.25*(yl(2)-yl(1));
yl_new = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];
ylim(yl_new)

% show asterisks for p values
for i = 1:size(overlap,2)
    text(i,mean(overlap(:,i))+0.02,get_asterisks(pval(i),1),'horizontalalignment','center',...
        'fontsize',20)
end

% asterisk for diff
plot([1 2],[ybar ybar],'k-','linewidth',2)
text(1.5,yast2,get_asterisks(pval_diff,1),'horizontalalignment','center',...
        'fontsize',20)

xticks([1 2])
xticklabels(xtext)
ylabel(ytext)
set(gca,'fontsize',15)

end