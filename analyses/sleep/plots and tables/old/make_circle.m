function make_circle

%alpha = randn(60,1)*.4+pi/2;
%alpha = pi/2;
% note that histogram shows thing ENDING at that time, not starting
times = [0 6*3600 12*3600 18*3600]';
polar = convert_times_to_polar(times);
alpha = polar;
figure
subplot(2,2,1)
circ_plot(alpha,'pretty','ro',true,'linewidth',2,'color','r'),
title('pretty plot style')
subplot(2,2,2)
circ_plot(alpha,'hist',[],20,true,true,'linewidth',2,'color','r')
title('hist plot style')
subplot(2,2,3)
circ_plot(alpha,[],'s')
title('non-fancy plot style')

%{
N = 256;                        % Nr Colours
c = colormap(jet(N));           % Define ‘colormap’
%th = linspace(0, 2*pi, N);      % Create Polar Grid

th = ones(1,N);
th(1:100) = 0.5;

r = 0.5:0.5:1;
[TH,R] = meshgrid(th,r);
[X,Y] = pol2cart(TH,R);
C = bsxfun(@times,(X + Y),th);  % 

figure(1)                       % Plot Flat Spiral Without Edges
hs = surf(X,Y,C);
view([0  90])
axis square
grid off
set(hs, 'EdgeColor','none')
%}

end