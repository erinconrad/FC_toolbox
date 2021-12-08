function make_error

thing1 = (1:10)';
thing2 = (1:10)';

figure
h = tiledlayout(1,3,'tilespacing','tight','padding','tight');

nexttile
stackedplot([thing1,thing2]);

ax2 = nexttile;
plot(thing1);
 % if I comment this out the error goes away

nexttile
stackedplot([thing1,thing2]);
title('test 3') 

xlabel(ax2,'x')
%{
Canvas update iteration limit exceeded. This can occur
if the scene is marked dirty during a drawnow.

Error in
matlab.graphics.chart.StackedLineChart/set.Title

Error in matlab.graphics.chart.Chart/title

Error in title (line 53)
            title(ax,args{:});

Error in make_error (line 18)
title('test 3')
%}



end