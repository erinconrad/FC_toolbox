function reproduce_problem
figure
set(gcf,'position',[1 1 500 500])
t = tiledlayout(1,2,"TileSpacing",'tight','padding','tight');

nexttile
plot(1,1,'o')

tt = tiledlayout(t,1,1,'tilespacing','none','padding','none');
tt.Layout.Tile = 2;
tt.Layout.TileSpan = [1 1];
ax1 = axes(tt);
plot(ax1,1,1,'o')
ax1.XAxisLocation = 'top';
ax1.YAxisLocation = 'left';
ax1.XTick = 1:10;
ax1.XTickLabel = arrayfun(@(x) sprintf('hello test test test %d',x),1:10,'UniformOutput',false);
ax1.XLim = [1 10];

ax2 = axes(tt);
ax2.XAxisLocation = 'bottom';
ax2.YAxisLocation = 'right';
plot(ax2,1,1,'o')
ax2.XTick = 1:10;
ax2.XTickLabel = arrayfun(@(x) sprintf('hello test test test %d',x),1:10,'UniformOutput',false);
ax2.XLim = [1 10];
ax2.Color = 'none';
ax2.Box = 'off';

end