function [mp,stp] = shaded_error_bars_fc(times,m,st,color)

if isempty(color)
    color = [0,0.4470, 0.7410];
end

%% Flip dimensions if needed
if size(times,1) > 1
    times = times';
end

if size(m,1) > 1
    m = m';
end

if size(st,1) >2
    st = st';
end

%% Plot the line
mp = plot(times,m,'color',color,'linewidth',3);
hold on

%% Plot the patch
upper = st(2,:);
lower = st(1,:);
in_between = [upper, fliplr(lower)];
x2 = [times,fliplr(times)];
x2(x2==inf) = nan;

nan_idx = isnan(in_between) | isnan(x2);

stp = fill(x2(~nan_idx), in_between(~nan_idx),color,'linestyle','none');
alpha(stp,0.4);

end