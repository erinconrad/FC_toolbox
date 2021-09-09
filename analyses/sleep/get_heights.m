function heights = get_heights(yl,pairs_to_plot)

if ~isempty(pairs_to_plot)
    levels = height_sub_analysis(pairs_to_plot);
    levels = [levels; levels(end)+1;levels(end)+2];
else
    levels = [1;2];
end
nlevels = length(unique(levels)); % +1 to add a bar for the main comparison

level_heights1 = yl(1) + (yl(2)-yl(1))*(1+0.1*(1:nlevels)');
level_heights2 = yl(1) + (yl(2)-yl(1))*(1.04+0.1*(1:nlevels)');

heights = nan(length(levels),2);

for h = 1:size(heights,1)
    heights(h,:) = [level_heights1(levels(h)) level_heights2(levels(h))];
end


end