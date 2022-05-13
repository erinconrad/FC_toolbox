function stats_out = paired_plot(data,ytext,xlabels,skip_legend,legend_loc)

if exist('skip_legend') == 0
    skip_legend = 0;
end

if exist('legend_loc') == 0
    legend_loc = [];
end

%% Parameters
pcolor = [0, 0.4470, 0.7410];
ncolor = [0.6350, 0.0780, 0.1840];
ecolor = [0.9290, 0.6940, 0.1250];

%% Define positive and negative
pos_diff = data(:,2) > data(:,1);
neg_diff = data(:,1) > data(:,2);
equal_diff = data(:,1) == data(:,2);

%% Sign rank test
[pval,~,stats] = signrank(data(:,1),data(:,2));
Tpos =stats.signedrank; % Tpos = positive-rank sum = sum of positive ranks

%% Do plot
pp = plot(data(pos_diff,1),data(pos_diff,2),'o','markeredgecolor',pcolor,'linewidth',2,...
    'MarkerFaceColor',pcolor);
hold on
np = plot(data(neg_diff,1),data(neg_diff,2),'^','markeredgecolor',...
    ncolor,'MarkerFaceColor',ncolor,'linewidth',2);
ep = plot(data(equal_diff,1),data(equal_diff,2),'s','markeredgecolor',...
    ecolor,'markerfacecolor',ecolor,'linewidth',2);

xlabel(sprintf('%s %s',ytext,sprintf(xlabels{1})))
ylabel(sprintf('%s %s',ytext,sprintf(xlabels{2})))


all_min = min([ylim,xlim]);
all_max = max([xlim,ylim]);
plot([all_min all_max],[all_min all_max],'k--','linewidth',2)


xlim([all_min all_max])
ylim([all_min all_max])
xl = xlim;
yl = ylim;

px = xl(1) + 0.01*(xl(2)-xl(1));
py = yl(1) + 0.99*(yl(2)-yl(1));
text(px,py,get_p_text(pval),'verticalalignment','top','fontsize',15)
legtext1 = sprintf('Higher %s',xlabels{2});
legtext2 = sprintf('Lower %s',xlabels{2});
legtext3 = 'Equal';

if ~skip_legend
    if isempty(legend_loc)
        legend([pp;np;ep],{legtext1,legtext2,legtext3},...
            'location','southeast','fontsize',15)
    else
        legend([pp;np;ep],{legtext1,legtext2,legtext3},...
            'position',legend_loc,'fontsize',15)
    end

end
set(gca,'fontsize',15)

stats_out.medians = [nanmedian(data(:,1)) nanmedian(data(:,2))];
stats_out.iqrs = [iqr(data(:,1))  iqr(data(:,2))];
stats_out.Tpos = Tpos;
stats_out.pval = pval;
stats_out.nhigher_n = [sum(data(:,1) < data(:,2)) size(data,1)];


end