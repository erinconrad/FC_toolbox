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
        errorbar(1:size(data,1),nanmedian(data,2),...
            nanmedian(data,2)-pr(:,1),pr(:,2)-nanmedian(data,2),'o','linewidth',2);
    case 'all'
        for j = 1:size(data,1)
            plot(j+0.05*randn(size(data,2),1),...
                data(j,:)','o','linewidth',2)
            hold on
        end
    case 'lines'
        for j = 1:size(data,2)
            plot([1 2],[data(1,j) data(2,j)],'k')
            hold on
                
        end
end
hold on
xlim([0 size(data,1)+1])
xticks(1:size(data,1))
xticklabels(xlabels)
xtickangle(30);
ylabel(ytext)
yl = ylim;
[p,post_hoc_p,which_groups]=get_and_plot_non_para_stats(yl,data,p_or_unp);
set(gca,'fontsize',15)

end