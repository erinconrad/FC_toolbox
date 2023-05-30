function stats = unpaired_plot(a,b,xlabs,ylab,which_test)

na = length(a);
nb = length(b);
plot(1+randn(na,1)*0.05,a,'o','linewidth',2)
hold on
plot(2+randn(nb,1)*0.05,b,'o','linewidth',2)
xticks([1 2])
xticklabels(xlabs)
ylabel(ylab)

switch which_test
    case 'non_para'
        [p,~,stats_stuff] = ranksum(a,b);
        stats.p = p;
        stats.ns = [sum(~isnan(a)) sum(~isnan(b))];
        stats.ranksum = stats_stuff.ranksum;
    case 'para'
        [~,p,~,ostats] = ttest2(a,b);
        stats.p = p;
        stats.tstat = ostats.tstat;
        stats.df = ostats.df;
        stats.means = [nanmean(a) nanmean(b)];
        stats.sd = [nanstd(a) nanstd(b)];
end
yl = ylim;
ybar = yl(1)+(yl(2)-yl(1))*1.1;
ytext = yl(1)+(yl(2)-yl(1))*1.15;
ylnew = [yl(1) yl(1)+(yl(2)-yl(1))*1.2];
ylim(ylnew)
plot([1 2],[ybar,ybar],'k','linewidth',2)
text(1.5,ytext,get_p_text(p),'fontsize',15,'horizontalalignment','center')
set(gca,'fontsize',15)

end