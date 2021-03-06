function [pm,c] = pretty_matrix(thing,labels,cutoff,clabel,oned)

pm = turn_nans_gray(thing);
xticklabels([])
yticklabels([])
hold on
plot(xlim,[cutoff cutoff],'k','linewidth',5)
if ~oned
    plot([cutoff cutoff],ylim,'k','linewidth',5)
end

xl = xlim;
new_xl = [xl(1) - 0.12*(xl(2)-xl(1)),xl(2)];
line_pos = xl(1) - 0.06*(xl(2)-xl(1));

top_text_y = (1+cutoff)/2;
bottom_text_y = (cutoff+size(thing,1))/2;

text(line_pos,top_text_y,labels{1},'horizontalalignment','center',...
    'rotation',90,'fontsize',20)%,'verticalalignment','middle')

text(line_pos,bottom_text_y,labels{2},'horizontalalignment','center',...
    'rotation',90,'fontsize',20)%,'verticalalignment','middle')
xlim(new_xl)


if ~oned
    c = colorbar('location','westoutside');
    ylabel(c,clabel,'fontsize',20)
end
set(gca,'fontsize',20)


end