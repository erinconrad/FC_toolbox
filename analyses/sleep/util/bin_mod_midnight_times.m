function [which_bins,nbins,edges] = bin_mod_midnight_times(mod_midnight)

bin_size = 60*10; % 10 minute bins
midnight = 24*3600;

edges = 0:bin_size:midnight;
nbins = length(edges)-1;

[~,~,which_bins] = histcounts(mod_midnight,edges);

end