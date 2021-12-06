function make_error

thing1 = (1:10)';
thing2 = (1:10)';

figure
tiledlayout(3,2,'tilespacing','tight','padding','tight')

nexttile
stackedplot([thing1,thing2]);

nexttile
plot(thing1);
xlabel('x') % if I comment this out the error goes away

nexttile
stackedplot([thing1,thing2]);
title('test 3') % error occurs here
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