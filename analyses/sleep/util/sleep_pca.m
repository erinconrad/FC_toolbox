function sleep_pca(spike_bins)

npts = size(spike_bins,1);
nbins = size(spike_bins,2);

% Subtract mean
spike_bins = spike_bins - nanmean(spike_bins,2);

[coeff,score,latent] = pca(spike_bins,'Rows','complete');

figure
plot(coeff(:,1))
hold on
plot(coeff(:,2))
plot(coeff(:,3))

end