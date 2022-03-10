function stats_out = plot_paired_data(varargin)

data = varargin{1};
xlabels = varargin{2}(1:2);
if length(varargin{2})>2
    legtext = varargin{2}{3};
else
    legtext = 'in sleep';
end
ytext = varargin{3};
p_or_unp = varargin{4};

if length(varargin) >= 5
    plot_type = varargin{5};
else
    plot_type = 'errorbar';
end

if length(varargin) >= 6
    is_ranking = strcmp('ranking',varargin{6});
else
    is_ranking = 0;
end

if length(varargin) >= 7
    do_plot = varargin{7};
else
    do_plot = 1;
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
        %pcolor = [0.4660, 0.6740, 0.1880];
        pcolor = [0, 0.4470, 0.7410];
        ncolor = [0.6350, 0.0780, 0.1840];
        %ecolor = [0, 0.4470, 0.7410];
        %ecolor = [0.4660, 0.6740, 0.1880]; 
        ecolor = [0.9290, 0.6940, 0.1250];
        if is_ranking
            pos_diff = data(1,:) > data(2,:);
            neg_diff = data(2,:) > data(1,:);
        else
            pos_diff = data(2,:) > data(1,:);
            neg_diff = data(1,:) > data(2,:);
        end
        [pval,~,stats] = signrank(data(1,:)',data(2,:)');
        Tpos =stats.signedrank; % Tpos = positive-rank sum = sum of positive ranks
        if do_plot
            equal_diff = data(1,:) == data(2,:);
            pp = plot(data(1,pos_diff),data(2,pos_diff),'o','markeredgecolor',pcolor,'linewidth',2,...
                'MarkerFaceColor',pcolor);
            hold on
            np = plot(data(1,neg_diff),data(2,neg_diff),'^','markeredgecolor',...
                ncolor,'MarkerFaceColor',ncolor,'linewidth',2);
            ep = plot(data(1,equal_diff),data(2,equal_diff),'s','markeredgecolor',...
                ecolor,'markerfacecolor',ecolor,'linewidth',2);

            xlabel(sprintf('%s %s',ytext,sprintf(xlabels{1})))
            ylabel(sprintf('%s %s',ytext,sprintf(xlabels{2})))

            pause(0.3)

            all_min = min([ylim,xlim]);
            all_max = max([xlim,ylim]);



            plot([all_min all_max],[all_min all_max],'k--','linewidth',2)


            xlim([all_min all_max])
            ylim([all_min all_max])
            xl = xlim;
            yl = ylim;

            px = xl(1) + 0.01*(xl(2)-xl(1));
            py = yl(1) + 0.99*(yl(2)-yl(1));
            %ylim
            text(px,py,get_p_text(pval),'verticalalignment','top','fontsize',15)
            %ylim
            legtext1 = sprintf('Higher %s',legtext);
            legtext2 = sprintf('Lower %s',legtext);
            legtext3 = 'Equal';

            legend([pp;np;ep],{legtext1,legtext2,legtext3},...
                'location','southeast','fontsize',15,'box','off')
        end
        
        stats_out.medians = [nanmedian(data(1,:)) nanmedian(data(2,:))];
        stats_out.iqrs = [iqr(data(1,:))  iqr(data(2,:))];
        stats_out.Tpos = Tpos;
        stats_out.pval = pval;
        stats_out.nhigher_n = [sum(data(1,:) < data(2,:)) size(data,2)];
        
        
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

if do_plot
set(gca,'fontsize',15)
end



end