function [which_bins,nbins,edges] = bin_mod_midnight_times(mod_midnight,edges)

bin_size = 60*10; % 10 minute bins
midnight = 24*3600;

if isempty(edges)
    edges = 0:bin_size:midnight;
end
nbins = length(edges)-1;

[~,~,which_bins] = histcounts(mod_midnight,edges);

end