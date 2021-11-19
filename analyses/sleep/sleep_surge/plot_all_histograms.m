function plot_all_histograms(hists,times,names)

midpoint = size(hists,2)/2;

figure
tiledlayout(10,10,'tilespacing','tight','padding','tight')
for i = 1:size(hists,1)
    nexttile
    plot(times,hists(i,:))
    hold on
    plot([0 0],ylim,'--')
    title(names{i})
end

end