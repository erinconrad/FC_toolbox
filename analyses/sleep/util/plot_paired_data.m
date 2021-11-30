function plot_paired_data(varargin)

data = varargin{1};
xlabels = varargin{2};
ytext = varargin{3};
p_or_unp = varargin{4};

if length(varargin) == 5
    plot_type = varargin{5};
else
    plot_type = 'errorbar';
end

%(data,xlabels,ytext,p_or_unp,plot_type)

%errorbar(nanmean(data,2),nanstd(data,[],2),'o','linewidth',2);
pr = prctile(data,[25,75],2);

switch plot_type
    case 'errorbar'
        for j = 1:size(data,1)
            errorbar(j,nanmedian(data(j,:),2),...
                nanmedian(data(j,:),2)-pr(j,1),pr(j,2)-nanmedian(data(j,:),2),'o','linewidth',2);
            hold on
        end
    case 'all'
        for j = 1:size(data,1)
            plot(j+0.05*randn(size(data,2),1),...
                data(j,:)','o','linewidth',2)
            hold on
        end
    case 'lines'
        for j = 1:size(data,2)
            plot([1 2],[data(1,j) data(2,j)],'color',[0.5 0.5 0.5])
            hold on
                
        end
        plot([1 2],[nanmedian(data(1,:)) nanmedian(data(2,:))],'k','linewidth',2)
    case 'violin'
        violin_erin(data')
        hold on
    case 'scatter'
        pcolor = [0.4660, 0.6740, 0.1880];
        ncolor = [0.6350, 0.0780, 0.1840];
        ecolor = [0, 0.4470, 0.7410];
        pos_diff = data(2,:) > data(1,:);
        neg_diff = data(1,:) > data(2,:);
        equal_diff = data(1,:) == data(2,:);
        pp = plot(data(1,pos_diff),data(2,pos_diff),'+','color',pcolor,'linewidth',2);
        hold on
        np = plot(data(1,neg_diff),data(2,neg_diff),'x','color',ncolor,'linewidth',2);
        ep = plot(data(1,equal_diff),data(2,equal_diff),'o','color',ecolor,'linewidth',2);
        
        xlabel(sprintf('%s in %s',ytext,sprintf(xlabels{1})))
        ylabel(sprintf('%s in %s',ytext,sprintf(xlabels{2})))
        pval = signrank(data(1,:)',data(2,:)');
        pause(0.3)
        xl = xlim;
        yl = ylim;
        all_min = min([ylim,xlim]);
        all_max = max([xlim,ylim]);
        pause(0.3)
        text(xl(1),yl(2),get_p_text(pval),'verticalalignment','top','fontsize',15)
        plot([all_min all_max],[all_min all_max],'k--','linewidth',2)
        legend([pp;np;ep],{'Higher in sleep','Lower in sleep','Equal wake and sleep'},...
            'location','southeast','fontsize',15)
        
        
        
end

if ~strcmp(plot_type,'scatter')
    xlim([0 size(data,1)+1])
    xticks(1:size(data,1))
    xticklabels(xlabels)
    xtickangle(30);
    ylabel(ytext)
    yl = ylim;
    [p,post_hoc_p,which_groups]=get_and_plot_non_para_stats(yl,data,p_or_unp);
end
set(gca,'fontsize',15)

end