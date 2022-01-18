function [mp,stp] = shaded_error_bars_xy(xm,ym,xl,yl,xu,yu,color)

if isempty(color)
    color = [0,0.4470, 0.7410];
end

%% Plot the line
mp = plot(xm,ym,'color',color,'linewidth',3);
hold on

%% Plot the patch
in_between = [yu, fliplr(yl)];
x2 = [xu,fliplr(xl)];
x2(x2==inf) = nan;

nan_idx = isnan(in_between) | isnan(x2);

stp = fill(x2(~nan_idx), in_between(~nan_idx),color,'linestyle','none');
alpha(stp,0.4);

end