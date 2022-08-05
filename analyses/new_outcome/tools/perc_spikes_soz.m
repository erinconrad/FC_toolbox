function prop = perc_spikes_soz(spikes,soz)

assert(all(~isnan(soz)))
nan_spikes = isnan(spikes);

sp = spikes(~nan_spikes);
sz = soz(~nan_spikes);

prop = sum(sp(sz==1))/sum(sp)*100;

end