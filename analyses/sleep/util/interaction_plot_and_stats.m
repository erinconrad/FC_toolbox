function interaction_plot_and_stats(varargin)

fc= [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250];

thing = varargin{1};
legend_text = varargin{2};
ytext = varargin{3};
xtext = varargin{4};
show_leg = varargin{5};

if length(varargin) == 6
    plot_type = varargin{6};
else
    plot_type = 'errorbar';
end

% loop groups
lp = nan(size(thing,1),1);
xticklocs = [];
xtextreap = {};
yl = [0 0];
for i = 1:size(thing,1)
    
    pr = squeeze(prctile(thing(i,:,:),[25,75],2));
    
    switch plot_type
        case 'errorbar'
    
            lp(i) = errorbar([i-0.2 i+0.2],squeeze(nanmedian(thing(i,:,:),2)),...
                squeeze(nanmedian(thing(i,:,:),2))-pr(:,1),...
                pr(:,2)-squeeze(nanmedian(thing(i,:,:),2)),'o','linewidth',2);
        case 'violin'
            violin_erin(squeeze(thing(i,:,:)),'x',[i-0.2 i+0.2],...
                'facecolor',fc(i,:))
    end
    
    xticklocs = [xticklocs,i-0.2,i+0.2];
    xtextreap = [xtextreap,xtext{1},xtext{2}];
    %}
    hold on
    curr_yl = ylim;
    if curr_yl(1) < yl(1), yl(1) = curr_yl(1); end
    if curr_yl(2) > yl(2), yl(2) = curr_yl(2); end    
end
ylim(yl)
xticks(xticklocs)
xticklabels(xtextreap)
xlim([0.3 size(thing,1)+0.7])
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

% show names
if show_leg
    legend(lp,legend_text,'location','northeastoutside','fontsize',15);
else
    yl = ylim;
    ytext = yl(1) + 1.05*(yl(2)-yl(1));
    newyl = [yl(1) yl(1) + 1.1*(yl(2)-yl(1))];
    for i = 1:length(legend_text)
        text(i,ytext,legend_text{i},'horizontalalignment','center','fontsize',15);
    end
end


ylim(newyl)

end