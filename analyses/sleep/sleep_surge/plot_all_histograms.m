function plot_all_histograms(hists)

midpoint = size(hists,2)/2;

figure
tiledlayout(10,10,'tilespacing','tight','padding','tight')
for i = 1:size(hists,1)
    nexttile
    plot(hists(i,:))
    hold on
    plot([midpoint midpoint],ylim,'--')
end

end