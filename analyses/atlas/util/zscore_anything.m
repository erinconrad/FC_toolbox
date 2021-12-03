function zscore = zscore_anything(thing,dim)

zscore = (thing-nanmean(thing,dim))./nanstd(thing,[],dim);

end