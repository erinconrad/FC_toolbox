function counts = sz_counts_bins(szs,nbins,times)

counts = zeros(nbins,1);
for z = 1:size(szs,1)
    counts(szs(z,2)) = counts(szs(z,2)) + 1;
end

figure
plot(times,counts,'linewidth',2)
hold on
plot([0 0],ylim,'k--','linewidth',2)
xlabel('Hours surrounding sleep onset')
ylabel('# Seizures (across all patients)')
xlim([times(1) times(end)])
set(gca,'fontsize',15)
title('Seizure timing surrounding sleep')

end